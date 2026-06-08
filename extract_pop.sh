#!/bin/bash
#SBATCH --job-name=extract_pop_genes
#SBATCH --output=logs/extract_pop_genes_%A_%a.out
#SBATCH --error=logs/extract_pop_genes_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-4

set -euo pipefail

VAR_DIR="variants"

POP_LIST=(
    ACB ASW ESN GWD LWK MSL YRI
    CDX CHB CHS JPT KHV
    CEU FIN GBR IBS TSI
    BEB GIH ITU PJL STU
)

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

START=$(( SLURM_ARRAY_TASK_ID * 5 ))
END=$(( START + 4 ))

echo "=== ARRAY TASK $SLURM_ARRAY_TASK_ID processing POP indices $START to $END ==="

for (( i=START; i<=END && i<${#POP_LIST[@]}; i++ )); do
    POP=${POP_LIST[$i]}

    if [[ "ACB ASW ESN GWD LWK MSL YRI" =~ $POP ]]; then SUPERPOP="AFR"; fi
    if [[ "CDX CHB CHS JPT KHV" =~ $POP ]]; then SUPERPOP="EAS"; fi
    if [[ "CEU FIN GBR IBS TSI" =~ $POP ]]; then SUPERPOP="EUR"; fi
    if [[ "BEB GIH ITU PJL STU" =~ $POP ]]; then SUPERPOP="SAS"; fi

    POP_DIR="${VAR_DIR}/${SUPERPOP}/${POP}"
    MERGED_VCF="${POP_DIR}/${POP}_merged.vcf.gz"

    if [[ ! -f "$MERGED_VCF" ]]; then
        echo "Skipping $POP (no merged VCF)"
        continue
    fi

    echo "Processing population: $POP (Superpop: $SUPERPOP)"

    mkdir -p "${POP_DIR}/warfarin" "${POP_DIR}/codeine"

    for GENE in "${!GENES[@]}"; do
        REGION="${GENES[$GENE]}"

        if [[ "$GENE" == "CYP2C9" || "$GENE" == "VKORC1" || "$GENE" == "CYP4F2" || "$GENE" == "GGCX" || "$GENE" == "CALU" ]]; then
            DRUG="warfarin"
        else
            DRUG="codeine"
        fi

        OUT_DIR="${POP_DIR}/${DRUG}/${GENE}"
        mkdir -p "$OUT_DIR"

        OUT_VCF="${OUT_DIR}/${POP}_${GENE}.vcf.gz"

        echo "  Extracting $GENE → $OUT_VCF"

        bcftools view -r "$REGION" -Oz -o "$OUT_VCF" "$MERGED_VCF"
        bcftools index -f "$OUT_VCF"
    done
done

echo "=== ARRAY TASK $SLURM_ARRAY_TASK_ID COMPLETE ==="
