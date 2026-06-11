#!/usr/bin/env bash

set -euo pipefail

METADATA="metadata/igsr_samples.tsv"
OUTDIR="metadata/superpops"

mkdir -p "$OUTDIR"

# Requested sample counts per superpopulation
AFR_TOTAL=50
EUR_TOTAL=25
EAS_TOTAL=12
SAS_TOTAL=13

# Header
HEADER=$(head -n 1 "$METADATA")

########################################
# Helper: sample for one superpopulation
########################################
sample_superpop() {
    local SUPER=$1
    local TOTAL=$2
    shift 2
    local POPS=("$@")

    local SUPERDIR="$OUTDIR/$SUPER"
    mkdir -p "$SUPERDIR"

    local NUM_POPS=${#POPS[@]}
    local BASE=$(( TOTAL / NUM_POPS ))
    local REM=$(( TOTAL % NUM_POPS ))

    echo "Processing $SUPER → total $TOTAL, $NUM_POPS populations"
    echo "  Base per population: $BASE, remainder: $REM"

    local COUNT_INDEX=0
    for POP in "${POPS[@]}"; do
        local EXTRA=0
        if [ "$COUNT_INDEX" -lt "$REM" ]; then
            EXTRA=1
        fi
        local COUNT=$(( BASE + EXTRA ))
        COUNT_INDEX=$(( COUNT_INDEX + 1 ))

        echo "    Population $POP → $COUNT samples"

        local OUTFILE="$SUPERDIR/${POP}.tsv"
        echo -e "$HEADER" > "$OUTFILE"

        awk -F'\t' -v pop="$POP" 'NR>1 && $4==pop' "$METADATA" \
            | shuf -n "$COUNT" >> "$OUTFILE"
    done
}

########################################
# Define populations per superpopulation
########################################

# AFR populations
AFR_POPS=(YRI LWK GWD MSL ESN ACB ASW)

# EUR populations
EUR_POPS=(CEU TSI FIN GBR IBS)

# EAS populations
EAS_POPS=(CHB JPT CHS CDX KHV)

# SAS populations
SAS_POPS=(GIH PJL BEB STU ITU)

########################################
# Run sampling for each superpopulation
########################################

sample_superpop "AFR" "$AFR_TOTAL" "${AFR_POPS[@]}"
sample_superpop "EUR" "$EUR_TOTAL" "${EUR_POPS[@]}"
sample_superpop "EAS" "$EAS_TOTAL" "${EAS_POPS[@]}"
sample_superpop "SAS" "$SAS_TOTAL" "${SAS_POPS[@]}"

echo "Done."
