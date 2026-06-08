#!/bin/bash

#SBATCH --job-name=multiqc_all
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err


set -euo pipefail

# Clean previous combined report
rm -rf reports/combined_multiqc
mkdir -p reports/combined_multiqc

# Run MultiQC across all fastp reports
multiqc reports -o reports/combined_multiqc -f

echo "MultiQC aggregation complete."
echo "Combined report available in reports/combined_multiqc"
