#!/bin/bash

#SBATCH --job-name=fastp_stu
#SBATCH --cpus-per-task=10
#SBATCH --mem=8G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=0-1%2   # 2 samples, run both concurrently

set -euo pipefail


# Build list of STU sample directories
SAMPLE_DIRS=( $(find data/SAS/STU -mindepth 1 -maxdepth 1 -type d | sort) )
SAMPLEDIR="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}"

SAMPLE=$(basename "$SAMPLEDIR")

REPORTDIR="reports/SAS/STU/$SAMPLE"
TRIMDIR="$SAMPLEDIR/trimmed"

# Clean previous outputs
rm -rf "$REPORTDIR" "$TRIMDIR"
mkdir -p "$REPORTDIR" "$TRIMDIR"

FASTQS=($(ls "$SAMPLEDIR"/*.fastq.gz 2>/dev/null || true))

echo "Processing STU sample: $SAMPLE..."

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
