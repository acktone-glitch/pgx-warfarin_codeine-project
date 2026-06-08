#!/bin/bash
#SBATCH --job-name=gatk_hc_AFR
#SBATCH --output=/etc/ace-data/home/rnagaba/reset/logs/hc_%A_%a.out
#SBATCH --error=/etc/ace-data/home/rnagaba/reset/logs/hc_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

# -------------------- CONFIG --------------------
POP="AFR"   # change if you run other populations

BASE="/etc/ace-data/home/rnagaba/reset"
DATA_DIR="${BASE}/data/${POP}"
OUT_DIR="${BASE}/results/gatk/${POP}"
GVCF_DIR="${OUT_DIR}/gvcfs"
TMP_DIR="${OUT_DIR}/tmp"
LOG_DIR="${BASE}/logs"
SAMPLE_LIST="${OUT_DIR}/samples_${POP}.txt"

REF="${BASE}/reference/hg38.fa"

# -------------------- PREP DIRS --------------------
mkdir -p "${GVCF_DIR}" "${TMP_DIR}" "${LOG_DIR}"

# -------------------- CHECK REFERENCE --------------------
if [[ ! -f "${REF}" ]]; then
    echo "[ERROR] Reference FASTA not found: ${REF}" >&2
    exit 1
fi

if [[ ! -f "${REF}.fai" ]]; then
    echo "[ERROR] FASTA index (.fai) missing for: ${REF}" >&2
    exit 1
fi

DICT="${REF%.fa}.dict"
if [[ ! -f "${DICT}" ]]; then
    echo "[ERROR] Dictionary (.dict) missing for: ${REF}" >&2
    exit 1
fi

echo "[INFO] Reference: ${REF}"

# -------------------- GATK BINARY --------------------
if [[ -n "${CONDA_BASE:-}" && -x "${CONDA_BASE}/envs/annotation/bin/gatk" ]]; then
    GATK="${CONDA_BASE}/envs/annotation/bin/gatk"
else
    GATK="gatk"
fi

echo "[INFO] GATK: $(${GATK} --version 2>&1 | head -1)"

# -------------------- SAMPLE LIST --------------------
if [[ ! -f "${SAMPLE_LIST}" ]]; then
    echo "[INFO] Sample list not found, creating: ${SAMPLE_LIST}"
    find "${DATA_DIR}" -name "*.bam" -not -name "*.bai" | sort > "${SAMPLE_LIST}"
fi

if [[ ! -s "${SAMPLE_LIST}" ]]; then
    echo "[ERROR] Sample list is empty: ${SAMPLE_LIST}" >&2
    exit 1
fi

# -------------------- THIS ARRAY TASK --------------------
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "[ERROR] SLURM_ARRAY_TASK_ID is not set. Run as an array job." >&2
    exit 1
fi

BAM=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLE_LIST}" || true)

if [[ -z "${BAM}" ]]; then
    echo "[ERROR] No BAM for task ${SLURM_ARRAY_TASK_ID} in ${SAMPLE_LIST}" >&2
    exit 1
fi

if [[ ! -f "${BAM}" ]]; then
    echo "[ERROR] BAM file does not exist: ${BAM}" >&2
    exit 1
fi

SAMPLE=$(basename "${BAM}" .bam)
GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"

echo "[INFO] Task:   ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] Sample: ${SAMPLE}"
echo "[INFO] BAM:    ${BAM}"
echo "[INFO] GVCF:   ${GVCF}"

# -------------------- SKIP IF DONE --------------------
if [[ -f "${GVCF}" && -f "${GVCF}.tbi" ]]; then
    echo "[INFO] Output already exists and is indexed, skipping."
    exit 0
fi

# -------------------- BAM INDEX --------------------
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
    echo "[INFO] Indexing BAM..."
    "${GATK}" BuildBamIndex -I "${BAM}"
fi

# -------------------- HAPLOTYPECALLER --------------------
echo "[INFO] Running HaplotypeCaller..."


${GATK} --java-options "-Xmx28g -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${TMP_DIR}" \
    HaplotypeCaller \
    -R "${REF}" \
    -I "${BAM}" \
    -O "${GVCF}" \
    -ERC GVCF \
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
    --native-pair-hmm-threads "${SLURM_CPUS_PER_TASK}" \
    --verbosity ERROR

echo "[INFO] Done: ${SAMPLE}"
