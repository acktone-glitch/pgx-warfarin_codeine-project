#!/bin/bash
#SBATCH --job-name=merge_all
#SBATCH --output=logs/merge_all_%A.out
#SBATCH --error=logs/merge_all_%A.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G


set -euo pipefail

# -------- DIRECTORIES --------
VAR_DIR="variants"

# Ensure logs directory exists
mkdir -p logs

echo "=== STARTING POPULATION AND SUPERPOP MERGES ==="

###############################################
# 1. MERGE SAMPLES → POPULATION VCFs
###############################################
echo "=== MERGING POPULATIONS ==="

for SUPERPOP in AFR EAS EUR SAS; do
    SUPERPOP_DIR="${VAR_DIR}/${SUPERPOP}"

    echo "Processing SUPERPOP: ${SUPERPOP}"

    for POP_DIR in "${SUPERPOP_DIR}"/*; do
        [[ -d "$POP_DIR" ]] || continue

        POP=$(basename "$POP_DIR")
        echo "  Merging POP: ${POP}"

        # Collect sample VCFs (exclude merged files)
        mapfile -t VCF_LIST < <(find "$POP_DIR" -maxdepth 1 -name "*.vcf.gz" ! -name "*merged*" | sort)

        if [[ ${#VCF_LIST[@]} -eq 0 ]]; then
            echo "    No sample VCFs found for ${POP}, skipping."
            continue
        fi

        POP_MERGED="${POP_DIR}/${POP}_merged.vcf.gz"

        echo "    Output: ${POP_MERGED}"
        echo "    Samples: ${#VCF_LIST[@]}"

        # Overwrite old merged files
        rm -f "$POP_MERGED" "$POP_MERGED.csi" 2>/dev/null || true

        bcftools merge -Oz -o "$POP_MERGED" "${VCF_LIST[@]}"
        bcftools index -f "$POP_MERGED"

        echo "    Completed POP merge: ${POP}"
    done
done


###############################################
# 2. MERGE POPULATIONS → SUPERPOPULATION VCFs
###############################################
echo "=== MERGING SUPERPOPULATIONS ==="

for SUPERPOP in AFR EAS EUR SAS; do
    SUPERPOP_DIR="${VAR_DIR}/${SUPERPOP}"

    echo "Merging SUPERPOP: ${SUPERPOP}"

    # Collect population-level merged VCFs
    mapfile -t POP_MERGED_LIST < <(find "$SUPERPOP_DIR" -name "*_merged.vcf.gz" | sort)

    if [[ ${#POP_MERGED_LIST[@]} -eq 0 ]]; then
        echo "  No population merged VCFs found for ${SUPERPOP}, skipping."
        continue
    fi

    SUPERPOP_MERGED="${SUPERPOP_DIR}/${SUPERPOP}_merged.vcf.gz"

    echo "  Output: ${SUPERPOP_MERGED}"
    echo "  Populations: ${#POP_MERGED_LIST[@]}"

    # Overwrite old merged files
    rm -f "$SUPERPOP_MERGED" "$SUPERPOP_MERGED.csi" 2>/dev/null || true

    bcftools merge -Oz -o "$SUPERPOP_MERGED" "${POP_MERGED_LIST[@]}"
    bcftools index -f "$SUPERPOP_MERGED"

    echo "  Completed SUPERPOP merge: ${SUPERPOP}"
done

echo "=== ALL MERGES COMPLETED SUCCESSFULLY ==="
