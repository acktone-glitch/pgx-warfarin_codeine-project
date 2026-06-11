#!/bin/bash
#SBATCH --job-name=vc_YRI
#SBATCH --output=logs/vc_YRI_%A_%a.out
#SBATCH --error=logs/vc_YRI_%A_%a.err
#SBATCH --array=0-7%5
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G


set -euo pipefail

# -------- CONFIG --------
REF="reference/hg38"
DATA_DIR="data/AFR/YRI"
VAR_DIR="variants/AFR/YRI"

mkdir -p "$VAR_DIR" logs

# -------- GET LIST OF BAM FILES FOR YRI ONLY --------
mapfile -t BAM_FILES < <(find "$DATA_DIR" -type f -path "*/aligned/*.bam" | sort)

BAM="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

if [[ -z "$BAM" ]]; then
    echo "No BAM found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Extract sample name
REL_PATH=${BAM#"$DATA_DIR"/}
SAMPLE=$(echo "$REL_PATH" | cut -d/ -f2)

SAMPLE_VCF="${VAR_DIR}/${SAMPLE}.vcf.gz"

echo "Processing YRI sample: $SAMPLE"
echo "  BAM: $BAM"
echo "  VCF OUT: $SAMPLE_VCF"

# -------- VARIANT CALLING --------
rm -f "$SAMPLE_VCF" "$SAMPLE_VCF.csi" 2>/dev/null || true

bcftools mpileup -Ou -f "$REF" "$BAM" \
    | bcftools call -mv -Oz -o "$SAMPLE_VCF"

bcftools index -f "$SAMPLE_VCF"

# -------- UPDATE SAMPLE LIST --------
POP_LIST="${VAR_DIR}/YRI_samples.txt"
grep -qx "$SAMPLE" "$POP_LIST" 2>/dev/null || echo "$SAMPLE" >> "$POP_LIST"

echo "YRI sample $SAMPLE completed."
