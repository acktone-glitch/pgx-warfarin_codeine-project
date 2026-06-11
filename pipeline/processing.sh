#!/bin/bash

#SBATCH --job-name=fastp_QC
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err

set -euo pipefail

# Activate fastp environment

set +u
source activate fastp
set -u


FAILED_LOG="failed_fastp.log"
: > "$FAILED_LOG"

# Loop through all sample directories
find data -mindepth 3 -maxdepth 3 -type d | while read -r SAMPLEDIR; do
    SUPERPOP=$(echo "$SAMPLEDIR" | cut -d/ -f2)
    POP=$(echo "$SAMPLEDIR" | cut -d/ -f3)
    SAMPLE=$(basename "$SAMPLEDIR")

    REPORTDIR="reports/$SUPERPOP/$POP/$SAMPLE"
    TRIMDIR="$SAMPLEDIR/trimmed"

    mkdir -p "$REPORTDIR" "$TRIMDIR"

    FASTQS=($(ls "$SAMPLEDIR"/*.fastq.gz 2>/dev/null || true))

    if [[ ${#FASTQS[@]} -eq 2 ]]; then
        # Paired-end
        fastp \
            -i "${FASTQS[0]}" \
            -I "${FASTQS[1]}" \
            -o "$TRIMDIR/${SAMPLE}_R1.trimmed.fastq.gz" \
            -O "$TRIMDIR/${SAMPLE}_R2.trimmed.fastq.gz" \
            -h "$REPORTDIR/${SAMPLE}_fastp.html" \
            -j "$REPORTDIR/${SAMPLE}_fastp.json" \
            --thread 4
    elif [[ ${#FASTQS[@]} -eq 1 ]]; then
        # Single-end
        fastp \
            -i "${FASTQS[0]}" \
            -o "$TRIMDIR/${SAMPLE}.trimmed.fastq.gz" \
            -h "$REPORTDIR/${SAMPLE}_fastp.html" \
            -j "$REPORTDIR/${SAMPLE}_fastp.json" \
            --thread 4
    else
        echo "No FASTQs found in $SAMPLEDIR"
        continue
    fi

    # Check if outputs exist and are non-empty
    if [[ ! -s "$REPORTDIR/${SAMPLE}_fastp.json" ]] || [[ ! -s "$REPORTDIR/${SAMPLE}_fastp.html" ]]; then
        echo "$SAMPLE ($SUPERPOP/$POP) → fastp report missing" >> "$FAILED_LOG"
    fi
    if [[ ! -s "$TRIMDIR"/*.trimmed.fastq.gz ]]; then
        echo "$SAMPLE ($SUPERPOP/$POP) → trimmed FASTQ missing" >> "$FAILED_LOG"
    fi
done

# Switch to MultiQC environment
set +u
source activate multiqc2
set -u


# Run MultiQC per superpopulation
for SUPERPOP in $(ls data); do
    mkdir -p reports/$SUPERPOP
    multiqc reports/$SUPERPOP -o reports/$SUPERPOP
done

echo "fastp trimming + MultiQC complete."
echo "Failed samples logged in $FAILED_LOG (if any)."
