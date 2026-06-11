#!/bin/bash
#SBATCH --job-name=align_auto
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --array=0-96%10

set -euo pipefail
mkdir -p logs

#############################################
# REFERENCE (classic BWA index)
#############################################

REF="reference/hg38"

#############################################
# FIND ALL SAMPLE DIRECTORIES WITH TRIMMED FASTQs
#############################################

# This finds directories like:
# data/SUPERPOP/POP/SAMPLE
# ONLY if they contain a "trimmed" folder
mapfile -t SAMPLE_DIRS < <(
    find data -mindepth 3 -maxdepth 3 -type d -exec test -d "{}/trimmed" \; -print | sort
)

# Stop the array if index is out of range
if [[ $SLURM_ARRAY_TASK_ID -ge ${#SAMPLE_DIRS[@]} ]]; then
    echo "Array index ${SLURM_ARRAY_TASK_ID} exceeds sample count ${#SAMPLE_DIRS[@]}"
    exit 0
fi

SAMPLEDIR="${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}"

SUPERPOP=$(echo "$SAMPLEDIR" | cut -d/ -f2)
POP=$(echo "$SAMPLEDIR" | cut -d/ -f3)
SAMPLE=$(basename "$SAMPLEDIR")

TRIMDIR="$SAMPLEDIR/trimmed"
ALIGNDIR="$SAMPLEDIR/aligned"

mkdir -p "$ALIGNDIR"

R1="${TRIMDIR}/${SAMPLE}_R1.trimmed.fastq.gz"
R2="${TRIMDIR}/${SAMPLE}_R2.trimmed.fastq.gz"

BAM="${ALIGNDIR}/${SAMPLE}.bam"

echo "----------------------------------------"
echo "Aligning sample: $SAMPLE"
echo "SUPERPOP: $SUPERPOP"
echo "POP:      $POP"
echo "R1:       $R1"
echo "R2:       $R2"
echo "Output:   $BAM"
echo "----------------------------------------"

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "ERROR: Missing trimmed FASTQs for $SAMPLE"
    exit 1
fi

#############################################
# ALIGN WITH CLASSIC BWA MEM
#############################################

bwa mem -t "$SLURM_CPUS_PER_TASK" "$REF" "$R1" "$R2" \
    | samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$BAM" -

samtools index "$BAM"

echo "Done: $SAMPLE"
