#!/bin/bash
#SBATCH --job-name=merge_superpops
#SBATCH --output=logs/merge_superpops_%A.out
#SBATCH --error=logs/merge_superpops_%A.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

###############################################
# CONFIG — CHANGE THIS TO YOUR PROJECT ROOT
###############################################
VAR_DIR="variants"

mkdir -p logs

echo "=== STARTING SUPERPOPULATION MERGES ==="

###############################################
# MERGE POPULATIONS → SUPERPOPULATION VCFs
###############################################

for SUPERPOP in AFR EAS EUR SAS; do
    SUPERPOP_DIR="${VAR_DIR}/${SUPERPOP}"

    echo ""
    echo "=== Processing SUPERPOP: ${SUPERPOP} ==="
    echo "Looking in: ${SUPERPOP_DIR}"

    # Collect population-level merged VCFs
    mapfile -t POP_MERGED_LIST < <(find "$SUPERPOP_DIR" -maxdepth 2 -name "*_merged.vcf.gz" | sort)

    echo "Found population merged VCFs:"
    printf '  %s\n' "${POP_MERGED_LIST[@]}"

    if [[ ${#POP_MERGED_LIST[@]} -eq 0 ]]; then
        echo "  No population merged VCFs found for ${SUPERPOP}, skipping."
        continue
    fi

    SUPERPOP_MERGED="${SUPERPOP_DIR}/${SUPERPOP}_merged.vcf.gz"

    echo "Output will be: ${SUPERPOP_MERGED}"
    echo "Number of populations: ${#POP_MERGED_LIST[@]}"

    # Overwrite old superpop merged files
    rm -f "$SUPERPOP_MERGED" "$SUPERPOP_MERGED.csi" 2>/dev/null || true

    # Merge populations
    bcftools merge -Oz -o "$SUPERPOP_MERGED" "${POP_MERGED_LIST[@]}"

    # Index
    bcftools index -f "$SUPERPOP_MERGED"

    echo "Completed SUPERPOP merge: ${SUPERPOP}"
done

echo ""
echo "=== ALL SUPERPOP MERGES COMPLETED SUCCESSFULLY ==="
