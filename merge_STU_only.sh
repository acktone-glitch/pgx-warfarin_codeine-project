#!/bin/bash
#SBATCH --job-name=merge_STU
#SBATCH --output=logs/merge_STU_%A.out
#SBATCH --error=logs/merge_STU_%A.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -euo pipefail

VAR_DIR="variants"
STU_DIR="${VAR_DIR}/SAS/STU"

mkdir -p logs

echo "=== MERGING STU POPULATION ==="

mapfile -t STU_VCFS < <(find "$STU_DIR" -maxdepth 1 -name "*.vcf.gz" ! -name "*merged*" | sort)

if [[ ${#STU_VCFS[@]} -eq 0 ]]; then
    echo "ERROR: No STU sample VCFs found."
    exit 1
fi

STU_MERGED="${STU_DIR}/STU_merged.vcf.gz"

echo "Output: $STU_MERGED"
rm -f "$STU_MERGED" "$STU_MERGED.csi" 2>/dev/null || true

bcftools merge --force-single -Oz -o "$STU_MERGED" "${STU_VCFS[@]}"
bcftools index -f "$STU_MERGED"

echo "=== STU MERGE COMPLETE ==="
