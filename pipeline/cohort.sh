#!/bin/bash
#SBATCH --job-name=vc_merge
#SBATCH --output=logs/vc_merge_%j.out
#SBATCH --error=logs/vc_merge_%j.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G

set -euo pipefail

VAR_DIR="variants"
COHORT_VCF="${VAR_DIR}/cohort.vcf.gz"

echo "Starting cohort merge and population splits..."

# -------- COLLECT ALL SAMPLE VCFs --------
mapfile -t SAMPLE_VCFS < <(find "$VAR_DIR" -type f -name "*.vcf.gz" \
    -not -name "cohort.vcf.gz" \
    -not -name "*_samples.txt" \
    -not -name "*_pop.vcf.gz" \
    -not -name "*_superpop.vcf.gz" | sort)

if [[ ${#SAMPLE_VCFS[@]} -eq 0 ]]; then
    echo "No sample VCFs found. Did Script 1 run?"
    exit 1
fi

VCF_LIST_FILE="$(mktemp)"
printf "%s\n" "${SAMPLE_VCFS[@]}" > "$VCF_LIST_FILE"

# -------- MERGE INTO COHORT --------
echo "Merging cohort VCF..."
rm -f "$COHORT_VCF" "$COHORT_VCF.tbi" 2>/dev/null || true

bcftools merge -Oz -o "$COHORT_VCF" --file-list "$VCF_LIST_FILE"
bcftools index -f "$COHORT_VCF"

rm -f "$VCF_LIST_FILE"

echo "Cohort VCF created: $COHORT_VCF"

# -------- SPLIT BY SUPERPOP --------
echo "Splitting by super-population..."

for SUPERPOP_DIR in "$VAR_DIR"/*/; do
    SUPERPOP=$(basename "$SUPERPOP_DIR")
    SUPERPOP_LIST="${SUPERPOP_DIR}/${SUPERPOP}_samples.txt"

    [[ -f "$SUPERPOP_LIST" ]] || continue

    SUPERPOP_VCF="${SUPERPOP_DIR}/${SUPERPOP}.vcf.gz"

    echo "  Creating $SUPERPOP VCF..."
    bcftools view -S "$SUPERPOP_LIST" "$COHORT_VCF" -Oz -o "$SUPERPOP_VCF"
    bcftools index -f "$SUPERPOP_VCF"
done

# -------- SPLIT BY POPULATION --------
echo "Splitting by population..."

for SUPERPOP_DIR in "$VAR_DIR"/*/; do
    for POP_DIR in "$SUPERPOP_DIR"*/; do
        POP=$(basename "$POP_DIR")
        POP_LIST="${POP_DIR}/${POP}_samples.txt"

        [[ -f "$POP_LIST" ]] || continue

        POP_VCF="${POP_DIR}/${POP}.vcf.gz"

        echo "  Creating $POP VCF..."
        bcftools view -S "$POP_LIST" "$COHORT_VCF" -Oz -o "$POP_VCF"
        bcftools index -f "$POP_VCF"
    done
done

echo "Merge and split complete."
