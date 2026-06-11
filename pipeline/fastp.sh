#!/bin/bash

#SBATCH --job-name=fastp_array
#SBATCH --cpus-per-task=64
#SBATCH --mem=48G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=0-98%10   # 99 samples, 10 at a time

set -euo pipefail

FASTP="conda run -n fastp fastp"

# Build sample list dynamically
SAMPLE_DIRS=( $(find data -mindepth 3 -maxdepth 3 -type d | sort) )
SAMPLEDIR="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}"

SUPERPOP=$(echo "$SAMPLEDIR" | cut -d/ -f2)
POP=$(echo "$SAMPLEDIR" | cut -d/ -f3)
SAMPLE=$(basename "$SAMPLEDIR")

REPORTDIR="reports/$SUPERPOP/$POP/$SAMPLE"
TRIMDIR="$SAMPLEDIR/trimmed"

# Clean previous outputs
rm -rf "$REPORTDIR" "$TRIMDIR"
mkdir -p "$REPORTDIR" "$TRIMDIR"

FASTQS=($(ls "$SAMPLEDIR"/*.fastq.gz 2>/dev/null || true))

echo "Processing $SAMPLE ($SUPERPOP/$POP)..."

if [[ ${#FASTQS[@]} -eq 2 ]]; then
    $FASTP \
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
