#!/bin/bash
#SBATCH --job-name=regen_cyp2d6
#SBATCH --output=logs/regen_cyp2d6_%A.out
#SBATCH --error=logs/regen_cyp2d6_%A.err
#SBATCH --cpus-per-task=5
#SBATCH --mem=15G

set -euo pipefail

# CYP2D6 coordinates on GRCh38
REGION="chr22:42125962-42131236"

# Root directory of your variants
ROOT_DIR=/etc/ace-data/home/rnagaba/reset/variants

# Populations to regenerate
POPS=(YRI MSL LWK GWD KHV JPT CHS CHB CDX IBS GBR CEU STU PJL ITU GIH BEB)

for POP in "${POPS[@]}"; do
    # Detect superpopulation folder automatically
    SUPERPOP=$(find $ROOT_DIR -type d -name $POP -printf '%h\n' | xargs basename)

    INPUT_VCF=$ROOT_DIR/$SUPERPOP/$POP/${POP}_merged.vcf.gz
    OUTPUT_DIR=$ROOT_DIR/$SUPERPOP/$POP/codeine/CYP2D6
    OUTPUT_VCF=$OUTPUT_DIR/${POP}_CYP2D6.vcf.gz

    if [ ! -s "$INPUT_VCF" ]; then
        echo "❌ Missing merged VCF for $POP: $INPUT_VCF"
        continue
    fi

    echo "=== Extracting CYP2D6 for $POP ==="
    bcftools view -r $REGION $INPUT_VCF -Oz -o $OUTPUT_VCF
    bcftools index $OUTPUT_VCF

    # Quick sanity check: ensure header line exists
    if ! zgrep -q "^#CHROM" $OUTPUT_VCF; then
        echo "⚠️  $OUTPUT_VCF missing #CHROM header"
    else
        echo "✅ $OUTPUT_VCF regenerated and indexed"
    fi
done

echo "All requested CYP2D6 VCFs regenerated."
