#!/bin/bash
#SBATCH --job-name=fastp_repl
#SBATCH --output=logs/fastp_repl_%A_%a.out
#SBATCH --error=logs/fastp_repl_%A_%a.err
#SBATCH --array=0-55%10
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

set -euo pipefail
mkdir -p logs

FASTP="conda run -n fastp fastp"

#############################################
# REPLACEMENT SAMPLE LIST (unique, cleaned)
#############################################

SAMPLES=(
HG01880
HG01882
HG01886
HG02013
HG01890
NA20274
NA19922
NA20342
NA20351
NA20355
HG02977
HG02981
HG03157
HG03139
HG03267
NA19020
NA19028
HG02813
HG02814
HG02840
HG02878
HG02982
NA18488
NA18498
NA18499
NA06989
HG00171
HG00173
HG01501
HG01503
HG01504
HG00097
HG00105
HG00106
HG00125
HG00145
NA20503
NA20511
NA20507
HG00956
NA18530
HG00407
NA18946
NA19001
HG01596
NA20864
NA21092
HG02651
HG02691
HG03723
HG03790
HG03646
HG03672
HG03603
HG03814
HG04144
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

REPORTDIR="reports/$SUPERPOP/$POP/$SAMPLE"
TRIMDIR="$SAMPLEDIR/trimmed"

# Clean previous outputs
rm -rf "$REPORTDIR" "$TRIMDIR"
mkdir -p "$REPORTDIR" "$TRIMDIR"

FASTQS=($(ls "$SAMPLEDIR"/*.fastq.gz 2>/dev/null || true))

echo "Processing $SAMPLE ($SUPERPOP/$POP)..."

#############################################
# RUN FASTP
#############################################

if [[ ${#FASTQS[@]} -eq 2 ]]; then
    $FASTP \
        -i "${FASTQS[0]}" \
        -I "${FASTQS[1]}" \
        -o "$TRIMDIR/${SAMPLE}_R1.trimmed.fastq.gz" \
        -O "$TRIMDIR/${SAMPLE}_R2.trimmed.fastq.gz" \
        -h "$REPORTDIR/${SAMPLE}_fastp.html" \
        -j "$REPORTDIR/${SAMPLE}_fastp.json" \
        --thread $SLURM_CPUS_PER_TASK
else
    echo "ERROR: $SAMPLE does not have exactly 2 FASTQ files" >&2
fi

echo "Finished $SAMPLE"
