#!/bin/bash
#SBATCH --job-name=vep_superpop
#SBATCH --output=logs/vep_superpop_%A.out
#SBATCH --error=logs/vep_superpop_%A.err
#SBATCH --cpus-per-task=15
#SBATCH --mem=62G


set -euo pipefail

CACHE_DIR=/etc/ace-data/home/rnagaba/reset/.vep/cache
PLUGIN_DIR=/etc/ace-data/home/rnagaba/reset/.vep/plugins
CLINVAR_VCF=$CACHE_DIR/clinvar/clinvar.vcf.gz
DBNSFP=$CACHE_DIR/dbNSFP/dbNSFP5.3.1a_grch38.gz

ROOT=/etc/ace-data/home/rnagaba/reset/variants

SUPERPOPS=(AFR EUR EAS SAS)

for SP in "${SUPERPOPS[@]}"; do
    INPUT_VCF=$ROOT/$SP/${SP}_merged.vcf.gz
    OUTPUT_VCF=$ROOT/$SP/${SP}_merged.vep.vcf.gz
    SUMMARY=$ROOT/$SP/${SP}_merged.summary.html
    LOG=$ROOT/$SP/${SP}_merged.vep.log

    if [ ! -s "$INPUT_VCF" ]; then
        echo "Missing merged VCF for $SP"
        continue
    fi

    echo "=== Annotating $SP merged VCF ==="

    if ! vep \
        --cache \
        --offline \
        --dir_cache $CACHE_DIR \
        --dir_plugins $PLUGIN_DIR \
        --species homo_sapiens \
        --assembly GRCh38 \
        --vcf \
        --input_file $INPUT_VCF \
        --output_file $OUTPUT_VCF \
        --force_overwrite \
        --plugin LoFtool \
        --plugin dbNSFP,$DBNSFP,ALL \
        --custom $CLINVAR_VCF,ClinVar,vcf,exact,0,CLNSIG \
        --fork 15 \
        --stats_file $SUMMARY \
        2>&1 | tee $LOG; then
            echo "❌ Error annotating $SP — see $LOG"
            continue
    fi

    echo "✅ Finished $SP → $OUTPUT_VCF"
done

echo "All superpopulation merged VCFs annotated."
