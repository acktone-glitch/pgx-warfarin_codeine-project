#!/bin/bash

#SBATCH --job-name=fastp_regen_reports
#SBATCH --cpus-per-task=5
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=0-12%4   # 13 samples, 4 concurrent

set -euo pipefail

# List of sample directories (Superpop/Pop/Sampleid)
SAMPLE_DIRS=(
    "AFR/ASW/NA19625"
    "AFR/ASW/NA20298"
    "AFR/ESN/HG03363"
    "AFR/LWK/NA19463"
    "AFR/YRI/NA18510"
    "AFR/YRI/NA18853"
    "AFR/YRI/NA19093"
    "EUR/FIN/HG00268"
    "EUR/IBS/HG01536"
    "EUR/GBR/ALL"
    "EAS/CDX/HG01811"
    "EAS/CHB/NA18614"
    "EAS/KHV/HG01848"
    "SAS/GIH/NA21106"
)

TARGET="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}"

SAMPLEDIR="data/$TARGET"
SUPERPOP=$(echo "$TARGET" | cut -d/ -f1)
POP=$(echo "$TARGET" | cut -d/ -f2)
SAMPLE=$(basename "$TARGET")

TRIMDIR="$SAMPLEDIR/trimmed"
REPORTDIR="reports/$SUPERPOP/$POP/$SAMPLE"

mkdir -p "$REPORTDIR"

R1="$TRIMDIR/${SAMPLE}_R1.trimmed.fastq.gz"
R2="$TRIMDIR/${SAMPLE}_R2.trimmed.fastq.gz"

echo "Regenerating report for $SUPERPOP/$POP/$SAMPLE..."

if [[ -f "$R1" && -f "$R2" ]]; then
    fastp \
        -i "$R1" \
        -I "$R2" \
        -h "$REPORTDIR/${SAMPLE}_fastp.html" \
        -j "$REPORTDIR/${SAMPLE}_fastp.json" \
        --disable_adapter_trimming \
        --disable_quality_filtering \
        --disable_length_filtering \
        --thread 2
else
    echo "ERROR: Trimmed FASTQs not found for $SAMPLE" >&2
fi
