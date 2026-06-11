#!/bin/bash
#SBATCH --job-name=vep_genes
#SBATCH --output=logs/vep_genes_%A.out
#SBATCH --error=logs/vep_genes_%A.err
#SBATCH --cpus-per-task=15
#SBATCH --mem=64G

set -u

# Absolute paths to cache and plugins
CACHE_DIR=/etc/ace-data/home/rnagaba/reset/.vep/cache
PLUGIN_DIR=/etc/ace-data/home/rnagaba/reset/.vep/plugins
CLINVAR_VCF=$CACHE_DIR/clinvar/clinvar.vcf.gz
DBNSFP=$CACHE_DIR/dbNSFP/dbNSFP5.3.1a_grch38.gz
LOFTOOL_SCORES=$PLUGIN_DIR/LoFtool_scores.txt

# Find all gene-level VCFs under drug folders (codeine, warfarin, etc.)
VCF_LIST=$(find /etc/ace-data/home/rnagaba/reset/variants -type f -path "*/codeine/CYP2D6/*.vcf.gz")

# Loop through each gene VCF
for INPUT_VCF in $VCF_LIST; do
    if [ ! -s "$INPUT_VCF" ]; then
        echo "Skipping empty or missing file: $INPUT_VCF"
        continue
    fi

    OUTPUT_VCF=${INPUT_VCF%.vcf.gz}.vep.vcf.gz
    LOG_FILE=${INPUT_VCF%.vcf.gz}.vep.log
    SUMMARY_HTML=${INPUT_VCF%.vcf.gz}.summary.html

    echo "=== Running VEP on $INPUT_VCF ==="

    if ! vep \
       --cache \
       --offline \
       --dir_cache $CACHE_DIR \
       --dir_plugins $PLUGIN_DIR \
       --species homo_sapiens \
       --assembly GRCh38 \
       --vcf \
       --compress_output bgzip \
       --input_file $INPUT_VCF \
       --output_file $OUTPUT_VCF \
       --force_overwrite \
       --plugin LoFtool,$LOFTOOL_SCORES \
       --plugin dbNSFP,$DBNSFP,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,CADD_phred,REVEL_score,MutationTaster_pred,gnomAD_exomes_AF,gnomAD_genomes_AF,LRT_score \
       --custom $CLINVAR_VCF,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
       --fork 15 \
       --everything \
       --stats_file $SUMMARY_HTML \
       2>&1 | tee $LOG_FILE; then
        echo "ERROR:failed processing $INPUT_VCF. See $LOG_FILE for details."
        continue
    fi

    echo "=== Finished $INPUT_VCF → $OUTPUT_VCF ==="
done

echo "All gene-level VEP runs complete."
