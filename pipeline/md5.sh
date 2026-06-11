#!/bin/bash

#SBATCH --job-name=MD5_Check
#SBATCH --cpus-per-task=15
#SBATCH --mem=24G
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err


set -euo pipefail

REF_MD5="igsr_md5.txt"
FAILED="failed_downloads.txt"
PASSED="passed_files.txt"

: > "$FAILED"
: > "$PASSED"

find data -type f -name "*.fastq.gz" | while read -r file; do
    fname=$(basename "$file")
    expected=$(grep -F "$fname" "$REF_MD5" | awk '{print $NF}')

    if [[ -z "$expected" ]]; then
        echo "$fname → no reference MD5 found" >> "$FAILED"
        continue
    fi

    if [[ ! -s "$file" ]]; then
        echo "$fname → empty file" >> "$FAILED"
        continue
    fi

    actual=$(md5sum "$file" | awk '{print $1}')

    if [[ "$actual" == "$expected" ]]; then
        echo "$fname → OK" >> "$PASSED"
    else
        url=$(grep -r "$fname" data/ --include="*_link" | head -n 1 | awk -F: '{print $2}' | awk '{print $1}')
        outdir=$(dirname "$file")
        echo "$expected $fname $url $outdir" >> "$FAILED"
    fi
done

echo "MD5 check complete."
echo "Passed files listed in $PASSED"
echo "Failed files listed in $FAILED"
