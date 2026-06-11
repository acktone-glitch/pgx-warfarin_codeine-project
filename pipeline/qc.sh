#!/bin/bash

#SBATCH --job-name=fastp_loop
#SBATCH --cpus-per-task=64
#SBATCH --mem=48G
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err

set -euo pipefail

FASTP="conda run -n fastp fastp"
MULTIQC="conda run -n multiqc2 multiqc"

# Loop through all sample directories
for SAMPLEDIR in $(find data -mindepth 3 -maxdepth 3 -type d); do
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
    elif [[ ${#FASTQS[@]} -eq 1 ]]; then
        $FASTP \
            -i "${FASTQS[0]}" \
            -o "$TRIMDIR/${SAMPLE}.trimmed.fastq.gz" \
            -h "$REPORTDIR/${SAMPLE}_fastp.html" \
            -j "$REPORTDIR/${SAMPLE}_fastp.json" \
            --thread 4
    else
        echo "No FASTQs found in $SAMPLEDIR"
    fi
done

# Run MultiQC per superpopulation
for SUPERPOP in $(ls data); do
    mkdir -p reports/$SUPERPOP
    $MULTIQC reports/$SUPERPOP -o reports/$SUPERPOP -f
done

echo "fastp trimming + MultiQC complete."
