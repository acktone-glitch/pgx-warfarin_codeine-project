#!/bin/bash
#SBATCH --job-name=fastp_three
#SBATCH --output=logs/fastp_three_%A_%a.out
#SBATCH --error=logs/fastp_three_%A_%a.err
#SBATCH --array=0-2%3
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

set -euo pipefail
mkdir -p logs


#############################################
# THREE SAMPLES YOU WANT TO PROCESS
#############################################

SAMPLES=(
NA20351
NA18870
HG03603
)

#############################################
# SELECT SAMPLE
#############################################

SID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Locate the sample directory under data/
SAMPLEDIR=$(find data -type d -name "$SID" | head -n 1)

if [[ -z "$SAMPLEDIR" ]]; then
    echo "ERROR: Could not find directory for $SID"
    exit 1
fi

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

#############################################
# RUN FASTP
#############################################

if [[ ${#FASTQS[@]} -eq 2 ]]; then
    fastp \
        -i "${FASTQS[0]}" \
        -I "${FASTQS[1]}" \
        -o "$TRIMDIR/${SAMPLE}_R1.trimmed.fastq.gz" \
        -O "$TRIMDIR/${SAMPLE}_R2.trimmed.fastq.gz" \
        -h "$REPORTDIR/${SAMPLE}_fastp.html" \
        -j "$REPORTDIR/${SAMPLE}_fastp.json" \
        --thread $SLURM_CPUS_PER_TASK
else
    echo "ERROR: $SAMPLE does not have exactly 2 FASTQ files" >&2
fi

echo "Finished $SAMPLE"
