#!/bin/bash
# =============================================================================
# rescue_afr_cyp2d6.sh
# Rescues AFR CYP2D6 variant 2 (QUAL=81, MQ=19) by relaxing MQ to >=18.
#
# DECISION LOG:
#   Variant 1: chr22:42129770 QUAL=3.2 MQ=6  → DROPPED (too low quality)
#   Variant 2: chr22:42130715 QUAL=81  MQ=19 → RESCUED (MQ relaxed to >=18)
#
#   Rationale: QUAL=81 is strong evidence the variant is real.
#   MQ=19 is one point below previous threshold. CYP2D6 pseudogene
#   region routinely produces low MQ — this is a known technical artefact
#   not a sign of a bad variant call. Flagged in analysis table.
# =============================================================================

#SBATCH --job-name=rescue_afr_cyp2d6
#SBATCH --output=logs/rescue_afr_cyp2d6_%A.out
#SBATCH --error=logs/rescue_afr_cyp2d6_%A.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail
mkdir -p logs

CACHE_DIR=/etc/ace-data/home/rnagaba/reset/.vep/cache
PLUGIN_DIR=/etc/ace-data/home/rnagaba/reset/.vep/plugins
CLINVAR_VCF=$CACHE_DIR/clinvar/clinvar.vcf.gz
DBNSFP=$CACHE_DIR/dbNSFP/dbNSFP5.3.1a_grch38.gz

INPUT_VCF=/etc/ace-data/home/rnagaba/reset/variants/AFR/codeine/CYP2D6/AFR_CYP2D6.vcf.gz
DIR=$(dirname "$INPUT_VCF")
VEP_OUT=$DIR/AFR_CYP2D6.vep.vcf.gz
ANALYSIS_OUT=$DIR/AFR_CYP2D6.analysis.vcf.gz
SUMMARY=$DIR/AFR_CYP2D6.summary.html
LOG=$DIR/AFR_CYP2D6.vep.log

# Relaxed filter — MQ>=18 to capture QUAL=81 variant
FILTER_EXPR='QUAL>=5 & INFO/DP>=2 & INFO/MQ>=18 & INFO/AC>=1'

echo "============================================"
echo "AFR CYP2D6 Rescue — MQ relaxed to >=18"
echo "Filter: $FILTER_EXPR"
echo "============================================"

echo "Raw variants before filtering:"
bcftools query -f '%CHROM\t%POS\tQUAL=%QUAL\tDP=%INFO/DP\tMQ=%INFO/MQ\tAC=%INFO/AC\n' "$INPUT_VCF"

echo ""
echo "STEP 1: VEP annotation..."
vep \
    --cache --offline \
    --dir_cache  "$CACHE_DIR" \
    --dir_plugins "$PLUGIN_DIR" \
    --species homo_sapiens \
    --assembly GRCh38 \
    --vcf \
    --input_file  "$INPUT_VCF" \
    --output_file "$VEP_OUT" \
    --force_overwrite \
    --plugin LoFtool \
    --plugin dbNSFP,"$DBNSFP",ALL \
    --custom "$CLINVAR_VCF",ClinVar,vcf,exact,0,CLNSIG \
    --fork 4 \
    --stats_file "$SUMMARY" \
    2>&1 | tee "$LOG"

echo "STEP 2: Filtering with $FILTER_EXPR ..."
bcftools filter \
    -i "$FILTER_EXPR" \
    "$VEP_OUT" \
| bcftools annotate \
    -x FORMAT/AD,FORMAT/PL \
    -Oz -o "$ANALYSIS_OUT" \
&& bcftools index "$ANALYSIS_OUT"

N_OUT=$(bcftools view -H "$ANALYSIS_OUT" | wc -l)
echo ""
echo "Variants rescued: $N_OUT"
echo "Rescued variant details:"
bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\tQUAL=%QUAL\tDP=%INFO/DP\tMQ=%INFO/MQ\tAC=%INFO/AC\n' \
    "$ANALYSIS_OUT"

echo ""
echo "Done: $ANALYSIS_OUT ready for transfer"
