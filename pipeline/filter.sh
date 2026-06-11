#!/bin/bash
# =============================================================================
# filter_vcfs.sh
# Filters all VEP-annotated VCFs for quality and strips unneeded fields.
# Run this on the HPC before transferring files to your laptop.
#
# What this does:
#   1. Applies quality filters (QUAL, DP, MQ, AC)
#   2. Removes FORMAT fields you don't need (AD, PL) to shrink file size
#   3. Outputs *.analysis.vcf.gz — one per gene per population
#
# Usage:
#   sbatch filter_vcfs.sh
#   OR run interactively:
#   bash filter_vcfs.sh
# =============================================================================

#SBATCH --job-name=filter_vcfs
#SBATCH --output=logs/filter_vcfs_%A.out
#SBATCH --error=logs/filter_vcfs_%A.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -e
mkdir -p logs

ROOT=/etc/ace-data/home/rnagaba/reset/variants

# ---------------------------------------------------------------------------
# Find all VEP-annotated filtered VCFs
# Expected path: ROOT/SUPERPOP/DRUG/GENE/SUPERPOP_GENE.filtered.vcf.gz
# ---------------------------------------------------------------------------
mapfile -t VCF_ARRAY < <(find "$ROOT" -maxdepth 4 -type f \
    \( -path "*/codeine/*/*.vep.vcf.gz" \
    -o -path "*/warfarin/*/*.vep.vcf.gz" \))

if [ ${#VCF_ARRAY[@]} -eq 0 ]; then
    echo "ERROR: No filtered VCFs found under $ROOT"
    exit 1
fi

echo "Found ${#VCF_ARRAY[@]} VCF(s) to filter"

PASS=0
FAIL=0

for VCF in "${VCF_ARRAY[@]}"; do

    echo "--------------------------------------------"
    echo "Input:  $VCF"

    # Build output path: replace .filtered.vcf.gz with .analysis.vcf.gz
    OUT="${VCF%.vcf.gz}.analysis.vcf.gz"
    echo "Output: $OUT"

    # -------------------------------------------------------------------------
    # FILTER EXPLAINED:
    #
    # -i 'QUAL>=10 & INFO/DP>=3 & INFO/MQ>=30 & INFO/AC>=1'
    #
    #   QUAL >= 10   → minimum confidence the variant is real
    #                  (below 10 is very likely a sequencing error)
    #
    #   INFO/DP >= 3 → at least 3 reads support this site
    #                  (1-2 reads is too low to trust in a population VCF)
    #
    #   INFO/MQ >= 30 → reads mapping with at least 99% confidence
    #                   (MQ < 30 means reads may be mapping to wrong location)
    #
    #   INFO/AC >= 1  → the ALT allele was actually observed in at least
    #                   one sample (removes sites where ALT was never seen)
    #
    # We do NOT filter on FILTER=="PASS" because your variants have
    # FILTER="." (no filter applied) — filtering on PASS would remove
    # everything.
    #
    # We do NOT apply strict AF cutoffs — rare population-specific
    # PGx variants matter and we want to keep them.
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    # FORMAT FIELD CLEANUP:
    #
    # -x FORMAT/AD  → removes per-allele read depths (not needed for pop PGx)
    # -x FORMAT/PL  → removes genotype likelihood scores (not needed)
    #
    # We KEEP FORMAT/GT (genotype calls) — essential for everything
    # -------------------------------------------------------------------------

    if bcftools filter \
        -i 'QUAL>=10 & INFO/DP>=3 & INFO/MQ>=30 & INFO/AC>=1' \
        "$VCF" \
    | bcftools annotate \
        -x FORMAT/AD,FORMAT/PL \
        -Oz -o "$OUT" \
    && bcftools index "$OUT"; then

        # Report how many variants passed
        N_IN=$(bcftools view -H "$VCF" | wc -l)
        N_OUT=$(bcftools view -H "$OUT" | wc -l)
        echo "Variants: $N_IN input → $N_OUT passed filter"
        echo "SUCCESS"
        ((PASS++)) || true
    else
        echo "FAILED"
        ((FAIL++)) || true
    fi

    echo ""
done

echo "============================================"
echo "Done: $PASS succeeded, $FAIL failed"
echo ""
echo "Files ready to transfer to laptop:"
find "$ROOT" -name "*.analysis.vcf.gz" | sort
echo "============================================"
