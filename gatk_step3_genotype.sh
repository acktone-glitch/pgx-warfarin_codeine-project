#!/bin/bash
#SBATCH --job-name=gatk_gt_AFR
#SBATCH --output=/etc/ace-data/home/rnagaba/reset/logs/gt_%j.out
#SBATCH --error=/etc/ace-data/home/rnagaba/reset/logs/gt_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# =============================================================================
#  STEP 3 — GenotypeGVCFs + hard filtering
#  Joint genotypes all samples, produces final filtered VCF
#
#  SUBMIT (only after step2 has COMPLETED successfully):
#    sbatch gatk_step3_genotype.sh
#
#  OUTPUT:
#    results/gatk/AFR/AFR_pgx_filtered.vcf.gz  ← use this for analysis
#
#  CHANGE POP= below for EAS / EUR / SAS
# =============================================================================

set -eo pipefail

POP="AFR"   # ← change this for other populations

BASE="/etc/ace-data/home/rnagaba/reset"
OUT_DIR="${BASE}/results/gatk/${POP}"
DB_DIR="${OUT_DIR}/genomicsdb"
LOG_DIR="${BASE}/logs"

mkdir -p "${OUT_DIR}/tmp" "${LOG_DIR}"

# ── REFERENCE ────────────────────────────────────────────────────────────────
REF="${BASE}/reference/hg38.fa"



echo "[INFO] Reference: ${REF}"

# Confirm .fai and .dict exist
if [[ ! -f "${REF}.fai" ]]; then
    echo "[ERROR] Missing ${REF}.fai — run: samtools faidx ${REF}"
    exit 1
fi
if [[ ! -f "${BASE}/reference/hg38.dict" ]]; then
    echo "[ERROR] Missing hg38.dict — run: gatk CreateSequenceDictionary -R ${REF}"
    exit 1
fi

# ── CONDA / GATK ─────────────────────────────────────────────────────────────
GATK="/etc/ace-data/home/rnagaba/.conda/envs/annotation/bin/gatk"
[[ ! -x "${GATK}" ]] && GATK="gatk"
echo "[INFO] GATK: $(${GATK} --version 2>&1 | head -1)"
# ── PATHS ────────────────────────────────────────────────────────────────────
DB_PATH="${DB_DIR}/pgx_${POP}"
RAW_VCF="${OUT_DIR}/${POP}_pgx_raw.vcf.gz"
SNP_RAW="${OUT_DIR}/${POP}_snps_raw.vcf.gz"
INDEL_RAW="${OUT_DIR}/${POP}_indels_raw.vcf.gz"
SNP_FILT="${OUT_DIR}/${POP}_snps_filtered.vcf.gz"
INDEL_FILT="${OUT_DIR}/${POP}_indels_filtered.vcf.gz"
FILTERED_VCF="${OUT_DIR}/${POP}_pgx_filtered.vcf.gz"

if [[ ! -d "${DB_PATH}" ]]; then
    echo "[ERROR] GenomicsDB not found at ${DB_PATH}. Did step2 complete?"
    exit 1
fi

# ── GENOTYPEGVCFS ─────────────────────────────────────────────────────────────
echo "[INFO] Running GenotypeGVCFs..."

${GATK} GenotypeGVCFs \
    --java-options "-Xmx28g -XX:ParallelGCThreads=4" \
    -R "${REF}" \
    -V "gendb://${DB_PATH}" \
    -O "${RAW_VCF}" \
    -L chr10:94762681-94855547 \
    -L chr10:96698415-96748853 \
    -L chr22:42126499-42130865 \
    -L chr19:15878663-15916641 \
    -L chr16:31096174-31111938 \
    -L chr2:85545546-85612953  \
    -L chr4:69048012-69118744  \
    -L chr6:154039662-154137408 \
    -L chr11:113409605-113475553 \
    -L chr7:128189484-128236610 \
    --tmp-dir "${OUT_DIR}/tmp" \
    --verbosity ERROR

echo "[INFO] Raw VCF written: ${RAW_VCF}"

# ── HARD FILTER SNPs ─────────────────────────────────────────────────────────
echo "[INFO] Filtering SNPs..."

${GATK} SelectVariants \
    --java-options "-Xmx16g" \
    -R "${REF}" -V "${RAW_VCF}" \
    --select-type-to-include SNP \
    -O "${SNP_RAW}"

${GATK} VariantFiltration \
    --java-options "-Xmx16g" \
    -R "${REF}" -V "${SNP_RAW}" \
    --filter-expression "QD < 2.0"              --filter-name "QD2" \
    --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
    --filter-expression "FS > 60.0"             --filter-name "FS60" \
    --filter-expression "SOR > 3.0"             --filter-name "SOR3" \
    --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O "${SNP_FILT}"

# ── HARD FILTER INDELs ───────────────────────────────────────────────────────
echo "[INFO] Filtering indels..."

${GATK} SelectVariants \
    --java-options "-Xmx16g" \
    -R "${REF}" -V "${RAW_VCF}" \
    --select-type-to-include INDEL \
    -O "${INDEL_RAW}"

${GATK} VariantFiltration \
    --java-options "-Xmx16g" \
    -R "${REF}" -V "${INDEL_RAW}" \
    --filter-expression "QD < 2.0"               --filter-name "QD2" \
    --filter-expression "FS > 200.0"             --filter-name "FS200" \
    --filter-expression "SOR > 10.0"             --filter-name "SOR10" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O "${INDEL_FILT}"

# ── MERGE ────────────────────────────────────────────────────────────────────
echo "[INFO] Merging SNPs + indels..."

${GATK} MergeVcfs \
    --java-options "-Xmx16g" \
    -I "${SNP_FILT}" \
    -I "${INDEL_FILT}" \
    -O "${FILTERED_VCF}"

# Clean up intermediates
rm -f "${SNP_RAW}" "${SNP_RAW}.tbi" \
      "${INDEL_RAW}" "${INDEL_RAW}.tbi" \
      "${SNP_FILT}" "${SNP_FILT}.tbi" \
      "${INDEL_FILT}" "${INDEL_FILT}.tbi"

# ── SUMMARY ──────────────────────────────────────────────────────────────────
echo ""
echo "=== DONE: ${POP} joint-called VCF ==="
echo "    Output: ${FILTERED_VCF}"
echo ""
echo "All variants:"
bcftools stats "${FILTERED_VCF}" | grep "^SN" | grep -E "records|SNPs|indels"
echo ""
echo "PASS only:"
bcftools view -f PASS "${FILTERED_VCF}" | bcftools stats | grep "^SN" | grep -E "records|SNPs|indels"
