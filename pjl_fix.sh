#!/bin/bash

#SBATCH --job-name=fastp_HG03619
#SBATCH --cpus-per-task=10
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err
#SBATCH --array=0-0%1   # single sample, run once

set -euo pipefail

# Build list with only HG03619 in SAS/PJL
SAMPLE_DIRS=( "data/SAS/PJL/HG03619" )
SAMPLEDIR="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}"

SAMPLE=$(basename "$SAMPLEDIR")

REPORTDIR="reports/SAS/PJL/$SAMPLE"
TRIMDIR="$SAMPLEDIR/trimmed"

# Clean previous outputs
rm -rf "$REPORTDIR" "$TRIMDIR"
mkdir -p "$REPORTDIR" "$TRIMDIR"

FASTQS=($(ls "$SAMPLEDIR"/*.fastq.gz 2>/dev/null || true))

echo "Processing PJL sample: $SAMPLE..."

if [[ ${#FASTQS[@]} -eq 2 ]]; then
    fastp \
        -i "${FASTQS[0]}" \
        -I "${FASTQS[1]}" \
        -o "$TRIMDIR/${SAMPLE}_R1.trimmed.fastq.gz" \
        -O "$TRIMDIR/${SAMPLE}_R2.trimmed.fastq.gz" \
        -h "$REPORTDIR/${SAMPLE}_fastp.html" \
        -j "$REPORTDIR/${SAMPLE}_fastp.json" \
        --thread 4
else
    echo "ERROR: $SAMPLE does not have exactly 2 FASTQ files" >&2
fi
