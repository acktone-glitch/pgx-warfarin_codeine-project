#!/bin/bash
#SBATCH --job-name=vc_array
#SBATCH --output=logs/vc_array_%A_%a.out
#SBATCH --error=logs/vc_array_%A_%a.err
#SBATCH --array=0-95%10
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G

set -euo pipefail

# -------- CONFIG --------
REF="reference/hg38"
DATA_DIR="data"
VAR_DIR="variants"

mkdir -p "$VAR_DIR" logs

# -------- GET LIST OF BAM FILES --------
mapfile -t BAM_FILES < <(find "$DATA_DIR" -type f -path "*/aligned/*.bam" | sort)

BAM="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

if [[ -z "$BAM" ]]; then
    echo "No BAM found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Example BAM path:
# data/AFR/GWD/HG02813/aligned/HG02813.bam
REL_PATH=${BAM#"$DATA_DIR"/}
SUPERPOP=$(echo "$REL_PATH" | cut -d/ -f1)
POP=$(echo "$REL_PATH" | cut -d/ -f2)
SAMPLE=$(echo "$REL_PATH" | cut -d/ -f3)

SAMPLE_DIR="${VAR_DIR}/${SUPERPOP}/${POP}"
mkdir -p "$SAMPLE_DIR"

SAMPLE_VCF="${SAMPLE_DIR}/${SAMPLE}.vcf.gz"

echo "Processing sample:"
echo "  SUPERPOP: $SUPERPOP"
echo "  POP:      $POP"
echo "  SAMPLE:   $SAMPLE"
echo "  BAM:      $BAM"
echo "  VCF OUT:  $SAMPLE_VCF"

# -------- VARIANT CALLING (overwrite old files) --------
rm -f "$SAMPLE_VCF" "$SAMPLE_VCF.tbi" 2>/dev/null || true

bcftools mpileup -Ou -f "$REF" "$BAM" \
    | bcftools call -mv -Oz -o "$SAMPLE_VCF"

bcftools index -f "$SAMPLE_VCF"

# -------- UPDATE SAMPLE LISTS --------
SUPERPOP_LIST="${VAR_DIR}/${SUPERPOP}/${SUPERPOP}_samples.txt"
POP_LIST="${VAR_DIR}/${SUPERPOP}/${POP}/${POP}_samples.txt"

mkdir -p "$(dirname "$SUPERPOP_LIST")"
mkdir -p "$(dirname "$POP_LIST")"

grep -qx "$SAMPLE" "$SUPERPOP_LIST" 2>/dev/null || echo "$SAMPLE" >> "$SUPERPOP_LIST"
grep -qx "$SAMPLE" "$POP_LIST" 2>/dev/null || echo "$SAMPLE" >> "$POP_LIST"

echo "Sample $SAMPLE completed."
