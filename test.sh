#!/bin/bash
#SBATCH --job-name=vep_test
#SBATCH --output=vep_test_%j.out
#SBATCH --error=vep_test_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

# Define paths (update these to your actual directories)
#CACHE_DIR=~/reset/.vep/cache
#PLUGIN_DIR=~/reset/.vep/plugins
#CLINVAR_VCF=$CACHE_DIR/clinvar/clinvar.vcf.gz
#DBNSFP=$PLUGIN_DIR/dbNSFP4.5a.gz-
#INPUT_VCF=variants/AFR/YRI/warfarin/CYP2C9/YRI_CYP2C9.vcf.gz
#OUTPUT_VCF=variants/AFR/YRI/warfarin/CYP2C9/YRI_CYP2C9.vep.vcf.gz

# Run VEP
vep \
   --cache \
   --offline \
   --dir_cache /etc/ace-data/home/rnagaba/reset/.vep/cache \
   --dir_plugins /etc/ace-data/home/rnagaba/reset/.vep/plugins \
   --species homo_sapiens \
   --assembly GRCh38 \
   --vcf \
   --input_file variants/AFR/YRI/warfarin/CYP2C9/YRI_CYP2C9.vcf.gz \
   --output_file variants/AFR/YRI/warfarin/CYP2C9/YRI_CYP2C9.vep.vcf.gz \
   --force_overwrite \
   --plugin LoFtool \
   --plugin dbNSFP,/etc/ace-data/home/rnagaba/reset/.vep/cache/dbNSFP/dbNSFP5.3.1a_grch38.gz,ALL \
   --custom /etc/ace-data/home/rnagaba/reset/.vep/cache/clinvar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG


echo "VEP test run complete. Output written to $OUTPUT_VCF"
