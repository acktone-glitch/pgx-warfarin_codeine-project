#!/bin/bash
#SBATCH --job-name=pharmcat_preprocess_raw
#SBATCH --output=logs/pharmcat_preprocess_raw_%A.out
#SBATCH --error=logs/pharmcat_preprocess_raw_%A.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

ROOT=/etc/ace-data/home/rnagaba/reset/variants
REF=/etc/ace-data/home/rnagaba/reset/reference/hg38.fa

SUPERPOPS=(AFR EUR EAS SAS)

# PharmCAT-supported genes (coordinates GRCh38)
declare -A GENES=(
  ["CYP2C9"]="chr10: 94938588-94990148"
  ["VKORC1"]="chr16: 31090842-31095980"
  ["CYP4F2"]="chr19: 15877893-15900528"
  ["GGCX"]="chr2: 85544720-85561819"
  ["CALU"]="chr7: 128739292-128773400"
  ["CYP2D6"]="chr22: 42125962-42131236"
#["OPRM1"]="chr6: 154010496-154246867"
  ["UGT2B7"]="chr4: 69051363-69112987"
  ["DRD2"]="chr11: 113409605-113476402"
)

for SP in "${SUPERPOPS[@]}"; do
    INPUT_VCF=$ROOT/$SP/${SP}_merged.vcf.gz
    OUTDIR=$ROOT/$SP/pharmcat_ready
    mkdir -p $OUTDIR

    echo "=== Processing $SP ==="

    # Extract only PharmCAT genes
    REGIONS=$(printf "%s " "${GENES[@]}")

    RAW_SUBSET=$OUTDIR/${SP}_pharmcat_subset.vcf.gz
    CLEAN_VCF=$OUTDIR/${SP}_pharmcat.vcf.gz

    echo "  → Extracting PharmCAT genes"
    bcftools view -r $REGIONS $INPUT_VCF -Oz -o $RAW_SUBSET
    bcftools index $RAW_SUBSET

    echo "  → Normalizing, splitting, removing annotations"
    bcftools annotate -x INFO,^FORMAT/GT $RAW_SUBSET \
      | bcftools norm -f $REF -m -both \
      | bcftools sort \
      | bgzip -c > $CLEAN_VCF

    bcftools index $CLEAN_VCF
    tabix -p vcf $CLEAN_VCF

    echo "  ✔️ PharmCAT-ready VCF created: $CLEAN_VCF"
done

echo "All superpopulations preprocessed for PharmCAT."
