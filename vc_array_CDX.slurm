#!/bin/bash
#SBATCH --job-name=vc_CDX
#SBATCH --output=logs/CDX_%A_%a.out
#SBATCH --error=logs/CDX_%A_%a.err
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G

set -euo pipefail

# -------- CONFIG --------
REF="reference/hg38"
DATA_DIR="data/EAS/CDX"
VAR_DIR="variants/EAS/CDX"

mkdir -p "$VAR_DIR" logs

# -------- GET LIST OF BAM FILES FOR CDX ONLY --------
mapfile -t BAM_FILES < <(find "$DATA_DIR" -type f -path "*/aligned/*.bam" | sort)

BAM="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

if [[ -z "$BAM" ]]; then
    echo "No BAM found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Extract sample name correctly
# data/EAS/CDX/HG00956/aligned/HG00956.bam
SAMPLE=$(basename "$BAM" .bam)

SAMPLE_VCF="${VAR_DIR}/${SAMPLE}.vcf.gz"

echo "Processing CDX sample: $SAMPLE"
echo "  BAM: $BAM"
echo "  VCF OUT: $SAMPLE_VCF"

# -------- VARIANT CALLING --------
rm -f "$SAMPLE_VCF" "$SAMPLE_VCF.csi" 2>/dev/null || true

bcftools mpileup -Ou -f "$REF" "$BAM" \
    | bcftools call -mv -Oz -o "$SAMPLE_VCF"

bcftools index -f "$SAMPLE_VCF"

# -------- UPDATE SAMPLE LIST --------
POP_LIST="${VAR_DIR}/CDX_samples.txt"
grep -qx "$SAMPLE" "$POP_LIST" 2>/dev/null || echo "$SAMPLE" >> "$POP_LIST"

echo "CDX sample $SAMPLE completed."
