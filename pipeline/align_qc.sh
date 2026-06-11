#!/bin/bash
#SBATCH --job-name=bam_qc
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-96%20   # adjust after counting BAMs

set -euo pipefail
mkdir -p logs reports

#############################################
# FIND ALL BAM FILES
#############################################

mapfile -t BAM_FILES < <(find data -type f -name "*.bam" | sort)

if [[ $SLURM_ARRAY_TASK_ID -ge ${#BAM_FILES[@]} ]]; then
    echo "Array index ${SLURM_ARRAY_TASK_ID} exceeds BAM count ${#BAM_FILES[@]}"
    exit 0
fi

BAM="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

# Extract SUPERPOP, POP, SAMPLE from path
SUPERPOP=$(echo "$BAM" | cut -d/ -f2)
POP=$(echo "$BAM" | cut -d/ -f3)
SAMPLE=$(basename "$BAM" .bam)

# Create hierarchical report directory
SAMPLEDIR="reports/${SUPERPOP}/${POP}/${SAMPLE}"
mkdir -p "$SAMPLEDIR"

echo "QC for $SAMPLE ($SUPERPOP/$POP)"

#############################################
# BASIC QC METRICS
#############################################

samtools flagstat "$BAM" > "${SAMPLEDIR}/${SAMPLE}.flagstat.txt"
samtools idxstats "$BAM" > "${SAMPLEDIR}/${SAMPLE}.idxstats.txt"
samtools stats "$BAM" > "${SAMPLEDIR}/${SAMPLE}.stats.txt"

echo "QC done for $SAMPLE → $SAMPLEDIR"
