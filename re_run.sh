#!/bin/bash

#SBATCH --job-name=fastp_rerun_all
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=0-16%4   # 17 samples total, 4 concurrent

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
    "EUR/GBR/HGXXXX"   # replace HGXXXX with actual GBR sample IDs
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

REPORTDIR="reports/$SUPERPOP/$POP/$SAMPLE"
TRIMDIR="$SAMPLEDIR/trimmed"

# Clean previous outputs
rm -rf "$REPORTDIR" "$TRIMDIR"
mkdir -p "$REPORTDIR" "$TRIMDIR"

R1=$(ls $SAMPLEDIR/*_1.fastq.gz)
R2=$(ls $SAMPLEDIR/*_2.fastq.gz)

echo "Re-running fastp for $SUPERPOP/$POP/$SAMPLE..."

fastp \
    -i "$R1" \
    -I "$R2" \
    -o "$TRIMDIR/${SAMPLE}_R1.trimmed.fastq.gz" \
    -O "$TRIMDIR/${SAMPLE}_R2.trimmed.fastq.gz" \
    -h "$REPORTDIR/${SAMPLE}_fastp.html" \
    -j "$REPORTDIR/${SAMPLE}_fastp.json" \
    --thread 8

