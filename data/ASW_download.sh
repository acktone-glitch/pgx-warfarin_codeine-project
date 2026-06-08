#!/bin/bash


set -euo pipefail

LINKFILE="data/AFR/ASW_link"
BASEDIR="data/AFR/ASW"

mkdir -p "$BASEDIR"

current_sample=""

while IFS= read -r line; do
    # Skip empty lines
    [[ -z "$line" ]] && continue

    if [[ "$line" =~ ^HG[0-9]+ || "$line" =~ ^NA[0-9]+ ]]; then
        # This is a sample ID
        current_sample="$line"
        echo "Processing sample $current_sample"
        mkdir -p "$BASEDIR/$current_sample"
    elif [[ "$line" =~ ^ftp ]]; then
        # This is a FASTQ link
        if [[ -z "$current_sample" ]]; then
            echo "Error: found link before sample ID → $line"
            exit 1
        fi
        echo "  Downloading $line"
        wget -c -P "$BASEDIR/$current_sample" "$line"
    else
        echo "Skipping unrecognized line: $line"
    fi
done < "$LINKFILE"

echo "All downloads for ASW complete."
