#!/bin/bash
#SBATCH --job-name=gatk_vs_bcftools
#SBATCH --output=/etc/ace-data/home/rnagaba/reset/logs/compare_%j.out
#SBATCH --error=/etc/ace-data/home/rnagaba/reset/logs/compare_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -eo pipefail

BASE="/etc/ace-data/home/rnagaba/reset"
OUT="${BASE}/results/comparison"
LOG_DIR="${BASE}/logs"
mkdir -p "${OUT}" "${LOG_DIR}"

BCFTOOLS="/etc/ace-data/home/rnagaba/.conda/envs/annotation/bin/bcftools"

# PGx intervals — restrict bcftools VCF to same regions as GATK
REGIONS="chr10:94762681-94855547,chr10:96698415-96748853,chr22:42126499-42130865,chr19:15878663-15916641,chr16:31096174-31111938,chr2:85545546-85612953,chr4:69048012-69118744,chr6:154039662-154137408,chr11:113409605-113475553,chr7:128189484-128236610"

POPS=(AFR EAS EUR SAS)

for POP in "${POPS[@]}"; do

    echo ""
    echo "========================================"
    echo "  Comparing: ${POP}"
    echo "========================================"

    GATK_VCF="${BASE}/results/gatk/${POP}/${POP}_pgx_filtered.vcf.gz"
    BCFT_VCF="${BASE}/variants/${POP}/${POP}_merged.vep.vcf.gz"
    OUT_POP="${OUT}/${POP}"
    mkdir -p "${OUT_POP}"

    # Check inputs
    if [[ ! -f "${GATK_VCF}" ]]; then
        echo "  WARN: GATK VCF missing — ${GATK_VCF}"
        # Try alternate path
        GATK_VCF="${BASE}/results/gatk/${POP}/${POP}_pgx_raw.vcf.gz"
        [[ ! -f "${GATK_VCF}" ]] && echo "  ERROR: No GATK VCF found for ${POP}" && continue
    fi
    if [[ ! -f "${BCFT_VCF}" ]]; then
        echo "  WARN: bcftools VCF not found at ${BCFT_VCF}"
        BCFT_VCF="${BASE}/variants/${POP}/${POP}_merged.vcf.gz"
        [[ ! -f "${BCFT_VCF}" ]] && echo "  ERROR: No bcftools VCF found for ${POP}" && continue
    fi

    echo "  GATK:      ${GATK_VCF}"
    echo "  bcftools:  ${BCFT_VCF}"

    # Index if needed
    for VCF in "${GATK_VCF}" "${BCFT_VCF}"; do
        if [[ ! -f "${VCF}.tbi" && ! -f "${VCF}.csi" ]]; then
            echo "  Indexing $(basename ${VCF})..."
            ${BCFTOOLS} index --tbi "${VCF}" || ${BCFTOOLS} index --csi "${VCF}"
        fi
    done

    # Step 1: Counts (GATK = PASS only, bcftools = PGx regions only)
    echo "  [1] Variant counts..."
    GATK_N=$(${BCFTOOLS} view -f PASS "${GATK_VCF}" 2>/dev/null | \
             ${BCFTOOLS} stats | grep "^SN" | grep "number of records" | awk '{print $NF}')
    BCFT_N=$(${BCFTOOLS} view -r "${REGIONS}" "${BCFT_VCF}" 2>/dev/null | \
             ${BCFTOOLS} stats | grep "^SN" | grep "number of records" | awk '{print $NF}')
    echo "    GATK PASS variants (PGx regions):      ${GATK_N}"
    echo "    bcftools variants  (PGx regions):      ${BCFT_N}"

    # Step 2: Extract bcftools PGx region subset to tmp file
    echo "  [2] Extracting PGx region from bcftools VCF..."
    BCFT_PGX="${OUT_POP}/${POP}_bcftools_pgx.vcf.gz"
    ${BCFTOOLS} view -r "${REGIONS}" -O z -o "${BCFT_PGX}" "${BCFT_VCF}"
    ${BCFTOOLS} index --tbi "${BCFT_PGX}"

    # Step 3: isec — unique and shared sites
    echo "  [3] Running bcftools isec..."
    ISEC_DIR="${OUT_POP}/isec"
    rm -rf "${ISEC_DIR}"
    mkdir -p "${ISEC_DIR}"

    ${BCFTOOLS} isec \
        -p "${ISEC_DIR}" \
        -O z \
        "${BCFT_PGX}" \
        "${GATK_VCF}" \
        2>/dev/null || true

    # Count results
    UNIQUE_BCFT=0; UNIQUE_GATK=0; SHARED=0
    [[ -f "${ISEC_DIR}/0000.vcf.gz" ]] && \
        UNIQUE_BCFT=$(${BCFTOOLS} view "${ISEC_DIR}/0000.vcf.gz" | grep -vc "^#" || echo 0)
    [[ -f "${ISEC_DIR}/0001.vcf.gz" ]] && \
        UNIQUE_GATK=$(${BCFTOOLS} view "${ISEC_DIR}/0001.vcf.gz" | grep -vc "^#" || echo 0)
    [[ -f "${ISEC_DIR}/0002.vcf.gz" ]] && \
        SHARED=$(${BCFTOOLS} view "${ISEC_DIR}/0002.vcf.gz" | grep -vc "^#" || echo 0)

    echo "    Unique to bcftools:      ${UNIQUE_BCFT}"
    echo "    Unique to GATK:          ${UNIQUE_GATK}"
    echo "    Shared by both:          ${SHARED}"

    # Step 4: Ti/Tv
    echo "  [4] Ti/Tv ratios..."
    TITV_GATK=$(${BCFTOOLS} stats -f PASS "${GATK_VCF}" 2>/dev/null | \
                grep "^TSTV" | awk '{print $5}')
    TITV_BCFT=$(${BCFTOOLS} stats "${BCFT_PGX}" 2>/dev/null | \
                grep "^TSTV" | awk '{print $5}')
    echo "    GATK Ti/Tv:              ${TITV_GATK:-NA}"
    echo "    bcftools Ti/Tv:          ${TITV_BCFT:-NA}"
    echo "    (Expected coding ~2.8)"

    # Step 5: SNP / indel counts
    echo "  [5] SNP/indel counts..."
    for ITEM in "GATK:${GATK_VCF}" "bcftools:${BCFT_PGX}"; do
        LABEL="${ITEM%%:*}"; VCF="${ITEM#*:}"
        SNPS=$(${BCFTOOLS} stats "${VCF}" 2>/dev/null | grep "^SN" | grep "SNPs" | awk '{print $NF}')
        INDELS=$(${BCFTOOLS} stats "${VCF}" 2>/dev/null | grep "^SN" | grep "indels" | awk '{print $NF}')
        echo "    ${LABEL}: SNPs=${SNPS:-0}  indels=${INDELS:-0}"
    done

    # Step 6: Summary TSV
    SUMMARY="${OUT_POP}/${POP}_comparison_summary.tsv"
    {
      printf "Metric\tValue\n"
      printf "Population\t%s\n"             "${POP}"
      printf "GATK_PASS_variants\t%s\n"     "${GATK_N}"
      printf "bcftools_PGx_variants\t%s\n"  "${BCFT_N}"
      printf "Unique_to_bcftools\t%s\n"     "${UNIQUE_BCFT}"
      printf "Unique_to_GATK\t%s\n"         "${UNIQUE_GATK}"
      printf "Shared\t%s\n"                 "${SHARED}"
      printf "GATK_TiTv\t%s\n"             "${TITV_GATK:-NA}"
      printf "bcftools_TiTv\t%s\n"         "${TITV_BCFT:-NA}"
    } > "${SUMMARY}"

    # Clean up tmp
    rm -f "${BCFT_PGX}" "${BCFT_PGX}.tbi"

    echo "  ✅ ${POP} complete → ${SUMMARY}"
done

# Cross-population summary table
echo ""
echo "========================================"
echo "  CROSS-POPULATION SUMMARY"
echo "========================================"
printf "%-6s %8s %8s %10s %10s %8s %8s %8s\n" \
    "Pop" "GATK" "bcftools" "Uniq_GATK" "Uniq_bcft" "Shared" "GATK_TiTv" "bcft_TiTv"
for POP in "${POPS[@]}"; do
    S="${OUT}/${POP}/${POP}_comparison_summary.tsv"
    [[ ! -f "${S}" ]] && continue
    G=$(grep "GATK_PASS"     "${S}" | cut -f2)
    B=$(grep "bcftools_PGx"  "${S}" | cut -f2)
    UG=$(grep "Unique_to_GA" "${S}" | cut -f2)
    UB=$(grep "Unique_to_bc" "${S}" | cut -f2)
    SH=$(grep "^Shared"      "${S}" | cut -f2)
    GT=$(grep "GATK_TiTv"    "${S}" | cut -f2)
    BT=$(grep "bcftools_TiTv" "${S}" | cut -f2)
    printf "%-6s %8s %8s %10s %10s %8s %8s %8s\n" \
        "${POP}" "${G}" "${B}" "${UG}" "${UB}" "${SH}" "${GT}" "${BT}"
done

echo ""
echo "=== All populations complete. Results in ${OUT}/ ==="
