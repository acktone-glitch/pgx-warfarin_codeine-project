#!/bin/bash
#SBATCH --job-name=variants_pipeline
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

set -euo pipefail

# -------- CONFIG --------
REF="reference/hg38"          # path to your reference
DATA_DIR="data"
VAR_DIR="variants"

mkdir -p "$VAR_DIR" logs

echo "Starting variant pipeline..."
echo "Reference: $REF"
echo "Data dir: $DATA_DIR"
echo "Variants dir: $VAR_DIR"

# -------- 1. PER-SAMPLE VARIANT CALLING --------
echo "Step 1: per-sample variant calling"

# Find all BAMs in the expected structure
mapfile -t BAM_FILES < <(find "$DATA_DIR" -type f -path "*/aligned/*.bam" | sort)

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "No BAM files found under $DATA_DIR. Check your paths."
    exit 1
fi

ALL_SAMPLE_VCFS=()

for BAM in "${BAM_FILES[@]}"; do
    # Example path:
    # data/AFR/GWD/HG02813/aligned/HG02813.bam
    REL_PATH=${BAM#"$DATA_DIR"/}          # AFR/GWD/HG02813/aligned/HG02813.bam
    SUPERPOP=$(echo "$REL_PATH" | cut -d/ -f1)
    POP=$(echo "$REL_PATH" | cut -d/ -f2)
    SAMPLE=$(echo "$REL_PATH" | cut -d/ -f3)

    SAMPLE_VAR_DIR="${VAR_DIR}/${SUPERPOP}/${POP}"
    mkdir -p "$SAMPLE_VAR_DIR"

    SAMPLE_VCF="${SAMPLE_VAR_DIR}/${SAMPLE}.vcf.gz"

    echo "Calling variants for: SUPERPOP=$SUPERPOP POP=$POP SAMPLE=$SAMPLE"
    echo "  BAM: $BAM"
    echo "  VCF: $SAMPLE_VCF"

    # bcftools variant calling
    bcftools mpileup -Ou -f "$REF" "$BAM" \
        | bcftools call -mv -Oz -o "$SAMPLE_VCF"

    bcftools index -f "$SAMPLE_VCF"

    ALL_SAMPLE_VCFS+=("$SAMPLE_VCF")

    # -------- build sample lists --------
    # super-population list
    SUPERPOP_LIST="${VAR_DIR}/${SUPERPOP}/${SUPERPOP}_samples.txt"
    if ! grep -qx "$SAMPLE" "$SUPERPOP_LIST" 2>/dev/null; then
        echo "$SAMPLE" >> "$SUPERPOP_LIST"
    fi

    # population list
    POP_LIST="${VAR_DIR}/${SUPERPOP}/${POP}/${POP}_samples.txt"
    if ! grep -qx "$SAMPLE" "$POP_LIST" 2>/dev/null; then
        echo "$SAMPLE" >> "$POP_LIST"
    fi
done

echo "Per-sample calling complete."
echo "Total samples: ${#ALL_SAMPLE_VCFS[@]}"

# -------- 2. MERGE INTO COHORT VCF --------
echo "Step 2: merging into cohort VCF"

COHORT_VCF="${VAR_DIR}/cohort.vcf.gz"

# Create a file list for bcftools merge
VCF_LIST_FILE="$(mktemp)"
printf "%s\n" "${ALL_SAMPLE_VCFS[@]}" > "$VCF_LIST_FILE"

bcftools merge -Oz -o "$COHORT_VCF" --file-list "$VCF_LIST_FILE"
bcftools index -f "$COHORT_VCF"

rm -f "$VCF_LIST_FILE"

echo "Cohort VCF: $COHORT_VCF"

# -------- 3. SPLIT BY SUPER-POPULATION --------
echo "Step 3: splitting by super-population"

for SUPERPOP_DIR in "$VAR_DIR"/*/; do
    SUPERPOP=$(basename "$SUPERPOP_DIR")
    # skip non-superpop dirs (like cohort.vcf.gz etc.)
    [[ "$SUPERPOP" == "cohort.vcf.gz" ]] && continue
    [[ "$SUPERPOP" == "cohort.af.vcf.gz" ]] && continue

    SUPERPOP_LIST="${SUPERPOP_DIR}/${SUPERPOP}_samples.txt"
    [[ -f "$SUPERPOP_LIST" ]] || continue

    SUPERPOP_VCF="${SUPERPOP_DIR}/${SUPERPOP}.vcf.gz"

    echo "Creating super-pop VCF for $SUPERPOP"
    bcftools view -S "$SUPERPOP_LIST" "$COHORT_VCF" -Oz -o "$SUPERPOP_VCF"
    bcftools index -f "$SUPERPOP_VCF"
done

# -------- 4. SPLIT BY POPULATION --------
echo "Step 4: splitting by population"

for SUPERPOP_DIR in "$VAR_DIR"/*/; do
    SUPERPOP=$(basename "$SUPERPOP_DIR")
    [[ "$SUPERPOP" == "cohort.vcf.gz" ]] && continue
    [[ "$SUPERPOP" == "cohort.af.vcf.gz" ]] && continue

    # iterate over POP subdirs
    for POP_DIR in "$SUPERPOP_DIR"*/; do
        [[ -d "$POP_DIR" ]] || continue
        POP=$(basename "$POP_DIR")

        POP_LIST="${POP_DIR}/${POP}_samples.txt"
        [[ -f "$POP_LIST" ]] || continue

        POP_VCF="${POP_DIR}/${POP}.vcf.gz"

        echo "Creating population VCF for $SUPERPOP / $POP"
        bcftools view -S "$POP_LIST" "$COHORT_VCF" -Oz -o "$POP_VCF"
        bcftools index -f "$POP_VCF"
    done
done

echo "All done."
echo "Structure now in: $VAR_DIR"
