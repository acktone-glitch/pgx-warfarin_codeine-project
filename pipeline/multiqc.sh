#!/bin/bash
#SBATCH --job-name=multiqc_reports
#SBATCH --output=logs/multiqc_reports_%A.out
#SBATCH --error=logs/multiqc_reports_%A.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -euo pipefail
mkdir -p logs

# Directory where all fastp reports live
REPORTS_DIR="reports"

# Output directory for MultiQC summary
OUTDIR="multiqc"
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

echo "Running MultiQC on all fastp reports under $REPORTS_DIR"

multiqc "$REPORTS_DIR" -o "$OUTDIR" --force

echo "MultiQC complete. Output stored in $OUTDIR/"
