#!/bin/bash

#SBATCH --job-name=SAS_Down_Samples
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=0-4

set -euo pipefail

# SAS populations in array order
POPS=(GIH PJL BEB STU ITU)

# Select population for this array task
POP=${POPS[$SLURM_ARRAY_TASK_ID]}

echo "Starting downloads for population: $POP"

LINKFILE="data/SAS/${POP}_link"
BASEDIR="data/SAS/${POP}"

mkdir -p "$BASEDIR"

current_sample=""

while IFS= read -r line; do
    [[ -z "$line" ]] && continue

    # Detect sample ID
    if [[ "$line" =~ ^HG[0-9]+ || "$line" =~ ^NA[0-9]+ ]]; then
        current_sample="$line"
        echo "Processing sample $current_sample"
        mkdir -p "$BASEDIR/$current_sample"

    # Detect FASTQ link
    elif [[ "$line" =~ ^ftp ]]; then
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

echo "All downloads for $POP complete."
