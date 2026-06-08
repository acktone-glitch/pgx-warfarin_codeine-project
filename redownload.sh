#!/bin/bash

#SBATCH --job-name=Redownload_Failed
#SBATCH --cpus-per-task=64
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err

set -euo pipefail

# Define sample → linkfile → basedir mapping
declare -A LINKFILES
declare -A BASEDIRS

# ACB samples
LINKFILES["HG02339"]="data/AFR/ACB/HG02339_links.txt"
BASEDIRS["HG02339"]="data/AFR/ACB/HG02339"

LINKFILES["HG02511"]="data/AFR/ACB/HG02511_links.txt"
BASEDIRS["HG02511"]="data/AFR/ACB/HG02511"

# TSI samples
LINKFILES["NA20775"]="data/EUR/TSI/NA20775_links.txt"
BASEDIRS["NA20775"]="data/EUR/TSI/NA20775"

LINKFILES["NA20771"]="data/EUR/TSI/NA20771_links.txt"
BASEDIRS["NA20771"]="data/EUR/TSI/NA20771"

# Loop through all failed samples
for SAMPLE in HG02339 HG02511 NA20775 NA20771; do
    LINKFILE="${LINKFILES[$SAMPLE]}"
    OUTDIR="${BASEDIRS[$SAMPLE]}"

    echo "Processing sample $SAMPLE"
    mkdir -p "$OUTDIR"

    while IFS= read -r url; do
        [[ -z "$url" ]] && continue
        echo "  Downloading $url"
        wget -c -P "$OUTDIR" "$url"
    done < "$LINKFILE"

    echo "Sample $SAMPLE complete."
done

echo "All failed samples re-downloaded."
