#!/bin/bash
#SBATCH --job-name=gatk_hc_EUR
#SBATCH --output=/etc/ace-data/home/rnagaba/reset/logs/hc_%A_%a.out
#SBATCH --error=/etc/ace-data/home/rnagaba/reset/logs/hc_%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G

set -eo pipefail

POP="EUR"
BASE="/etc/ace-data/home/rnagaba/reset"
OUT_DIR="${BASE}/results/gatk/${POP}"
GVCF_DIR="${OUT_DIR}/gvcfs"
TMP_DIR="${OUT_DIR}/tmp"
LOG_DIR="${BASE}/logs"
SAMPLE_LIST="${OUT_DIR}/samples_${POP}.txt"
REF="${BASE}/reference/hg38.fa"

mkdir -p "${GVCF_DIR}" "${TMP_DIR}" "${LOG_DIR}"

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

# ── THIS SAMPLE ───────────────────────────────────────────────────────────────
BAM=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${SAMPLE_LIST}")
[[ -z "${BAM}" ]] && { echo "[ERROR] No BAM for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

SAMPLE=$(basename "${BAM}" .bam)
GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"
RG_BAM="${TMP_DIR}/${SAMPLE}.rg.bam"

echo "[INFO] Task:   ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] Sample: ${SAMPLE}"
echo "[INFO] BAM:    ${BAM}"
echo "[INFO] GVCF:   ${GVCF}"

# Skip if already done
if [[ -f "${GVCF}" && -f "${GVCF}.tbi" ]]; then
    echo "[INFO] Already complete — skipping"
    exit 0
fi

# ── STEP A: ADD READ GROUPS ───────────────────────────────────────────────────
echo "[INFO] Adding read groups..."

${GATK} AddOrReplaceReadGroups \
    --java-options "-Xmx16g" \
    -I "${BAM}" \
    -O "${RG_BAM}" \
    --RGID "${SAMPLE}" \
    --RGSM "${SAMPLE}" \
    --RGLB "${SAMPLE}" \
    --RGPL ILLUMINA \
    --RGPU "${SAMPLE}" \
    --TMP_DIR "${TMP_DIR}" \
    --VALIDATION_STRINGENCY SILENT \
    --CREATE_INDEX true

echo "[INFO] Read groups added"

# ── STEP B: HAPLOTYPECALLER ───────────────────────────────────────────────────
echo "[INFO] Running HaplotypeCaller..."

${GATK} HaplotypeCaller \
    --java-options "-Xmx28g -XX:ParallelGCThreads=4" \
    -R "${REF}" \
    -I "${RG_BAM}" \
    -O "${GVCF}" \
    -ERC GVCF \
    --sample-name "${SAMPLE}" \
    -L chr10:94762681-94855547  \
    -L chr10:96698415-96748853  \
    -L chr22:42126499-42130865  \
    -L chr19:15878663-15916641  \
    -L chr16:31096174-31111938  \
    -L chr2:85545546-85612953   \
    -L chr4:69048012-69118744   \
    -L chr6:154039662-154137408 \
    -L chr11:113409605-113475553 \
    -L chr7:128189484-128236610 \
    --tmp-dir "${TMP_DIR}" \
    --native-pair-hmm-threads "${SLURM_CPUS_PER_TASK}" \
    --verbosity ERROR

echo "[INFO] GVCF written: ${GVCF}"

# ── CLEANUP ───────────────────────────────────────────────────────────────────
rm -f "${RG_BAM}" "${RG_BAM%.bam}.bai"
echo "[INFO] Done: ${SAMPLE}"
