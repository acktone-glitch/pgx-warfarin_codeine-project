#!/bin/bash
#SBATCH --job-name=regen_reports
#SBATCH --output=logs/regen_reports_%A_%a.out
#SBATCH --error=logs/regen_reports_%A_%a.err
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

set -euo pipefail
mkdir -p logs

#############################################
# SAMPLES MISSING REPORTS
#############################################

SAMPLES=(
NA19028
HG03646
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

TRIMDIR="$SAMPLEDIR/trimmed"
REPORTDIR="reports/$SUPERPOP/$POP/$SAMPLE"

mkdir -p "$REPORTDIR"

R1="$TRIMDIR/${SAMPLE}_R1.trimmed.fastq.gz"
R2="$TRIMDIR/${SAMPLE}_R2.trimmed.fastq.gz"

echo "Regenerating fastp reports for $SAMPLE ($SUPERPOP/$POP)..."

#############################################
# RUN FASTP (REPORTS ONLY)
#############################################

fastp \
    -i "$R1" \
    -I "$R2" \
    -h "$REPORTDIR/${SAMPLE}_fastp.html" \
    -j "$REPORTDIR/${SAMPLE}_fastp.json" \
    --disable_adapter_trimming \
    --dont_overwrite \
    --thread $SLURM_CPUS_PER_TASK

echo "Done with $SAMPLE"
