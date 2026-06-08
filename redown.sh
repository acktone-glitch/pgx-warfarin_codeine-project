#!/bin/bash

#SBATCH --job-name=MD5_Redownload
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=0-$(($(wc -l < failed_downloads.txt)-1))

set -euo pipefail

REF_MD5="igsr_md5.txt"
FAILED_LIST="failed_downloads.txt"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$FAILED_LIST")

EXPECTED=$(echo "$LINE" | awk '{print $1}')
FNAME=$(echo "$LINE" | awk '{print $2}')
URL=$(echo "$LINE" | awk '{print $3}')
OUTDIR=$(echo "$LINE" | awk '{print $4}')

mkdir -p "$OUTDIR"

echo "Re-downloading $FNAME into $OUTDIR"
wget -c -O "$OUTDIR/$FNAME" "$URL"

ACTUAL=$(md5sum "$OUTDIR/$FNAME" | awk '{print $1}')
if [[ "$ACTUAL" == "$EXPECTED" ]]; then
    echo "$FNAME → OK after re-download"
else
    echo "$FNAME → STILL FAILED (expected $EXPECTED, got $ACTUAL)"
fi
