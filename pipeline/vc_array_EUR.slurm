#!/bin/bash
#SBATCH --job-name=vc_EUR
#SBATCH --output=logs/vc_EUR_%A_%a.out
#SBATCH --error=logs/vc_EUR_%A_%a.err
#SBATCH --array=0-24%5
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G

set -euo pipefail

REF="reference/hg38"
DATA_DIR="data/EUR"
VAR_DIR="variants/EUR"

mkdir -p "$VAR_DIR" logs

mapfile -t BAM_FILES < <(find "$DATA_DIR" -type f -path "*/aligned/*.bam" | sort)
BAM="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

if [[ -z "$BAM" ]]; then
    echo "No BAM found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

REL_PATH=${BAM#"$DATA_DIR"/}
POP=$(echo "$REL_PATH" | cut -d/ -f1)
SAMPLE=$(echo "$REL_PATH" | cut -d/ -f2)

SAMPLE_DIR="${VAR_DIR}/${POP}"
mkdir -p "$SAMPLE_DIR"

SAMPLE_VCF="${SAMPLE_DIR}/${SAMPLE}.vcf.gz"

rm -f "$SAMPLE_VCF" "$SAMPLE_VCF.csi" 2>/dev/null || true

bcftools mpileup -Ou -f "$REF" "$BAM" \
    | bcftools call -mv -Oz -o "$SAMPLE_VCF"

bcftools index -f "$SAMPLE_VCF"

POP_LIST="${SAMPLE_DIR}/${POP}_samples.txt"
grep -qx "$SAMPLE" "$POP_LIST" 2>/dev/null || echo "$SAMPLE" >> "$POP_LIST"
