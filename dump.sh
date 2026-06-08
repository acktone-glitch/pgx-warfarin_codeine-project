#!/bin/bash

# File: move_corrupted_samples.sh
# Purpose: Move corrupted sample directories into data/corrupt/ and metadata/corrupt/

set -euo pipefail

# Create corrupt directories if they don't exist
mkdir -p data/corrupt
mkdir -p metadata/corrupt

# List of corrupted samples (SUPERPOP/POP/SAMPLEID)
corrupted_samples=(
"AFR/ACB/HG02505"
"AFR/ACB/HG02339"
"AFR/ACB/HG02334"
"AFR/ACB/HG02323"
"AFR/ACB/HG02419"
"AFR/ASW/NA19713"
"AFR/ASW/NA19625"
"AFR/ASW/NA20298"
"AFR/ASW/NA20322"
"AFR/ASW/NA19707"
"AFR/ESN/HG03123"
"AFR/ESN/HG02973"
"AFR/ESN/HG03195"
"AFR/ESN/HG03135"
"AFR/ESN/HG03363"
"AFR/LWK/NA19311"
"AFR/LWK/NA19463"
"AFR/YRI/NA18510"
"AFR/YRI/NA18853"
"AFR/YRI/NA19093"
"AFR/YRI/NA19159"
"EUR/FIN/HG00268"
"EUR/FIN/HG00380"
"EUR/IBS/HG01710"
"EUR/IBS/HG01786"
"EUR/IBS/HG01536"
"EUR/GBR/HG00245"
"EUR/GBR/HG00099"
"EUR/GBR/HG00146"
"EUR/GBR/HG00122"
"EUR/GBR/HG00256"
"EUR/TSI/NA20775"
"EUR/TSI/NA20769"
"EUR/TSI/NA20771"
"EAS/CDX/HG01811"
"EAS/CHB/NA18614"
"EAS/CHS/HG00674"
"EAS/JPT/NA19075"
"EAS/JPT/NA18957"
"EAS/KHV/HG01848"
"SAS/GIH/NA21106"
"SAS/GIH/NA20871"
"SAS/PJL/HG03490"
"SAS/PJL/HG02601"
"SAS/ITU/HG03775"
"SAS/ITU/HG04118"
"SAS/STU/HG03856"
"SAS/STU/HG03894"
"SAS/BEB/HG03012"
"SAS/BEB/HG04183"
"SAS/BEB/HG03802"
)

echo "Moving corrupted sample directories..."

for entry in "${corrupted_samples[@]}"; do
    echo "Processing $entry"

    # Move data directory
    if [[ -d "data/$entry" ]]; then
        mkdir -p "data/corrupt/$(dirname "$entry")"
        mv "data/$entry" "data/corrupt/$entry"
    else
        echo "WARNING: data/$entry not found"
    fi

    # Move metadata directory
    if [[ -d "metadata/superpops/$entry" ]]; then
        mkdir -p "metadata/corrupt/$(dirname "$entry")"
        mv "metadata/superpops/$entry" "metadata/corrupt/$entry"
    else
        echo "WARNING: metadata/superpops/$entry not found"
    fi
done

echo "Done. All corrupted samples moved."

