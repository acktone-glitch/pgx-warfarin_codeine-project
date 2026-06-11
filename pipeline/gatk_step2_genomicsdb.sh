#!/bin/bash
#SBATCH --job-name=gatk_db_AFR
#SBATCH --output=/etc/ace-data/home/rnagaba/reset/logs/db_%j.out
#SBATCH --error=/etc/ace-data/home/rnagaba/reset/logs/db_%j.err

#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# =============================================================================
#  STEP 2 — GenomicsDBImport
#  Combines all per-sample GVCFs into a GenomicsDB workspace
#
#  SUBMIT (only after ALL step1 jobs have COMPLETED successfully):
#    sbatch gatk_step2_genomicsdb.sh
#
#  CHANGE POP= below for EAS / EUR / SAS
# =============================================================================

set -eo pipefail

POP="AFR"   # ← change this for other populations

BASE="/etc/ace-data/home/rnagaba/reset"
OUT_DIR="${BASE}/results/gatk/${POP}"
GVCF_DIR="${OUT_DIR}/gvcfs"
DB_DIR="${OUT_DIR}/genomicsdb"
LOG_DIR="${BASE}/logs"
SAMPLE_LIST="${OUT_DIR}/samples_${POP}.txt"

mkdir -p "${DB_DIR}" "${OUT_DIR}/tmp" "${LOG_DIR}"

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
# ── BUILD SAMPLE MAP ─────────────────────────────────────────────────────────
SAMPLE_MAP="${OUT_DIR}/sample_map_${POP}.tsv"
> "${SAMPLE_MAP}"

MISSING=0
while IFS= read -r BAM; do
    SAMPLE=$(basename "${BAM}" .bam)
    GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"
    if [[ -f "${GVCF}" && -f "${GVCF}.tbi" ]]; then
        echo -e "${SAMPLE}\t${GVCF}" >> "${SAMPLE_MAP}"
    else
        echo "[WARN] Missing GVCF for ${SAMPLE}"
        MISSING=$(( MISSING + 1 ))
    fi
done < "${SAMPLE_LIST}"

N=$(wc -l < "${SAMPLE_MAP}")
echo "[INFO] ${N} samples in map | ${MISSING} missing GVCFs"

if [[ "${N}" -eq 0 ]]; then
    echo "[ERROR] No GVCFs found. Did step1 complete?"
    exit 1
fi

# ── INTERVAL LIST ────────────────────────────────────────────────────────────
INTERVAL_FILE="${OUT_DIR}/pgx_intervals.list"
cat > "${INTERVAL_FILE}" << 'EOF'
chr10:94762681-94855547
chr10:96698415-96748853
chr22:42126499-42130865
chr19:15878663-15916641
chr16:31096174-31111938
chr2:85545546-85612953
chr4:69048012-69118744
chr6:154039662-154137408
chr11:113409605-113475553
chr7:128189484-128236610
EOF

# ── GENOMICSDBIMPORT ──────────────────────────────────────────────────────────
DB_PATH="${DB_DIR}/pgx_${POP}"

if [[ -d "${DB_PATH}" ]]; then
    echo "[WARN] Removing existing DB at ${DB_PATH}"
    rm -rf "${DB_PATH}"
fi

echo "[INFO] Building GenomicsDB at ${DB_PATH}..."

${GATK} GenomicsDBImport \
    --java-options "-Xmx28g -XX:ParallelGCThreads=4" \
    --sample-name-map "${SAMPLE_MAP}" \
    --genomicsdb-workspace-path "${DB_PATH}" \
    -L "${INTERVAL_FILE}" \
    --reader-threads "${SLURM_CPUS_PER_TASK}" \
    --batch-size 50 \
    --tmp-dir "${OUT_DIR}/tmp" \
    --verbosity ERROR

echo "[INFO] GenomicsDB complete: ${DB_PATH}"
