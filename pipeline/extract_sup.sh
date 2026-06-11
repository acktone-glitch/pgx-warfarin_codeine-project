#!/bin/bash
#SBATCH --job-name=extract_superpop_genes
#SBATCH --output=logs/extract_superpop_genes_%A.out
#SBATCH --error=logs/extract_superpop_genes_%A.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

VAR_DIR="variants"

declare -A GENES=(
    ["CYP2C9"]="chr10:94938488-94990248"
    ["VKORC1"]="chr16:31090742-31096080"
    ["CYP4F2"]="chr19:15877793-15900628"
    ["GGCX"]="chr2:85544620-85561919"
    ["CALU"]="chr7:128739192-128773500"
    ["CYP2D6"]="chr22:42125862-42131336"
    ["OPRM1"]="chr6:154010396-154246967"
    ["UGT2B7"]="chr4:69051263-69113087"
    ["DRD2"]="chr11:113409505-113476502"
)

SUPERPOPS=("AFR" "EAS" "EUR" "SAS")

echo "=== STARTING SUPERPOP-LEVEL EXTRACTION ==="

for SUPERPOP in "${SUPERPOPS[@]}"; do
    SUPERPOP_DIR="${VAR_DIR}/${SUPERPOP}"
    MERGED_VCF="${SUPERPOP_DIR}/${SUPERPOP}_merged.vcf.gz"

    if [[ ! -f "$MERGED_VCF" ]]; then
        echo "Skipping $SUPERPOP (no superpop merged VCF)"
        continue
    fi

    echo "Processing superpopulation: $SUPERPOP"

    mkdir -p "${SUPERPOP_DIR}/warfarin" "${SUPERPOP_DIR}/codeine"

    for GENE in "${!GENES[@]}"; do
        REGION="${GENES[$GENE]}"

        if [[ "$GENE" == "CYP2C9" || "$GENE" == "VKORC1" || "$GENE" == "CYP4F2" || "$GENE" == "GGCX" || "$GENE" == "CALU" ]]; then
            DRUG="warfarin"
        else
            DRUG="codeine"
        fi

        OUT_DIR="${SUPERPOP_DIR}/${DRUG}/${GENE}"
        mkdir -p "$OUT_DIR"

        OUT_VCF="${OUT_DIR}/${SUPERPOP}_${GENE}.vcf.gz"

        echo "  Extracting $GENE → $OUT_VCF"

        bcftools view -r "$REGION" -Oz -o "$OUT_VCF" "$MERGED_VCF"
        bcftools index -f "$OUT_VCF"
    done
done

echo "=== SUPERPOP-LEVEL EXTRACTION COMPLETE ==="
