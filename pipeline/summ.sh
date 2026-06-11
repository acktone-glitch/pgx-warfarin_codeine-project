#!/bin/bash
#SBATCH --job-name=vep_regen_summaries
#SBATCH --output=logs/vep_regen_summaries_%A.out
#SBATCH --error=logs/vep_regen_summaries_%A.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G


set -euo pipefail

CACHE_DIR=/etc/ace-data/home/rnagaba/reset/.vep/cache
PLUGIN_DIR=/etc/ace-data/home/rnagaba/reset/.vep/plugins

# List of CYP2D6 annotated VCFs missing summary.html
FILES=(
  "/etc/ace-data/home/rnagaba/reset/variants/AFR/GWD/codeine/CYP2D6/GWD_CYP2D6.vep.vcf.gz"
  "/etc/ace-data/home/rnagaba/reset/variants/AFR/LWK/codeine/CYP2D6/LWK_CYP2D6.vep.vcf.gz"
  "/etc/ace-data/home/rnagaba/reset/variants/AFR/MSL/codeine/CYP2D6/MSL_CYP2D6.vep.vcf.gz"
)

for FILE in "${FILES[@]}"; do
    if [ ! -s "$FILE" ]; then
        echo "Annotated file missing: $FILE"
        continue
    fi

    SUMMARY_HTML=${FILE%.vcf.gz}.summary.html
    LOG_FILE=${FILE%.vcf.gz}.regen.log

    echo "Regenerating summary for $FILE"

    vep \
      --cache \
      --offline \
      --dir_cache $CACHE_DIR \
      --dir_plugins $PLUGIN_DIR \
      --species homo_sapiens \
      --assembly GRCh38 \
      --vcf \
      --input_file $FILE \
      --output_file /dev/null \
      --stats_file $SUMMARY_HTML \
      --no_stats_update \
      2>&1 | tee $LOG_FILE

    echo "Summary regenerated: $SUMMARY_HTML"
done
