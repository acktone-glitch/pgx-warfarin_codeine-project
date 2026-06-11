#!/bin/bash
# =============================================================================
# parse_vep.sh
# Parses all VEP-annotated filtered VCFs for warfarin and codeine gene panels.
# Outputs per-VCF: variant table, CSQ annotation table, genotype matrix (CSV).
#
# What this script does, step by step:
#   1. Finds all *.filtered.vcf.gz files under variants/
#   2. For each VCF, launches an R session that:
#      a. Reads the VCF using VariantAnnotation
#      b. Extracts the variant-level INFO table
#      c. Parses the CSQ (VEP consequence) field into a clean table
#         - one row per VARIANT (takes the most severe transcript consequence)
#      d. Extracts the genotype matrix and cleans up sample names
#      e. Saves three separate CSVs: variants, csq, genotypes
# =============================================================================

#SBATCH --job-name=parse_vep
#SBATCH --output=logs/parse_vep_%A.out
#SBATCH --error=logs/parse_vep_%A.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

source ~/.bashrc
conda activate analysis

# ---------------------------------------------------------------------------
# STEP 1: Find all filtered VCFs
# We look 4 levels deep: ROOT/SUPERPOP/DRUG/GENE/file.filtered.vcf.gz
# ---------------------------------------------------------------------------
ROOT=/etc/ace-data/home/rnagaba/reset/variants

VCF_LIST=$(find "$ROOT" -maxdepth 4 -type f \
    \( -path "*/codeine/*/*.filtered.vcf.gz" -o -path "*/warfarin/*/*.filtered.vcf.gz" \))

mapfile -t VCF_ARRAY <<< "$VCF_LIST"

if [ ${#VCF_ARRAY[@]} -eq 0 ]; then
    echo "ERROR: No filtered VCFs found under $ROOT"
    echo "Expected paths like: $ROOT/AFR/warfarin/VKORC1/AFR_VKORC1.filtered.vcf.gz"
    exit 1
fi

echo "Found ${#VCF_ARRAY[@]} VCF(s) to process."

# ---------------------------------------------------------------------------
# STEP 2: Loop over each VCF and run the R parser
# ---------------------------------------------------------------------------
for VCF in "${VCF_ARRAY[@]}"; do

    echo ""
    echo "============================================"
    echo "Processing: $VCF"
    echo "============================================"

    export VCF_FILE="$VCF"

    Rscript - << 'REOF'

# =============================================================================
# R PARSING SCRIPT
# Reads a single VEP-annotated VCF and outputs three clean CSV files.
# =============================================================================

# ---------------------------------------------------------------------------
# Load libraries
# VariantAnnotation: Bioconductor package for reading VCF files
# stringr: for string manipulation (cleaning sample names, parsing CSQ)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(stringr)
})

# ---------------------------------------------------------------------------
# Get the VCF file path passed in from bash
# ---------------------------------------------------------------------------
vcf_file <- Sys.getenv("VCF_FILE")
if (vcf_file == "") stop("VCF_FILE environment variable not set.")

out_dir <- dirname(vcf_file)
# base name: e.g. "AFR_VKORC1" (strips path and .filtered.vcf.gz)
base    <- sub("\\.filtered\\.vcf\\.gz$", "", basename(vcf_file))

cat("Parsing:", vcf_file, "\n")
cat("Output directory:", out_dir, "\n")
cat("Base name:", base, "\n")

# ---------------------------------------------------------------------------
# STEP A: Read the VCF
#
# readVcf() loads the entire VCF into a VCF object with three parts:
#   - rowRanges: genomic coordinates (CHROM, POS, REF, ALT)
#   - info:      INFO field values (DP, AC, AN, CSQ, etc.)
#   - geno:      genotype data per sample (GT, AD, PL)
#
# NOTE: We do NOT filter on FILTER=="PASS" here.
# Your variants have FILTER="." which means "no filter applied" — this is
# normal for bcftools-called cohort VCFs. Filtering on PASS would wipe
# everything out. We keep all variants and you can filter later by
# consequence, AF, or quality scores.
# ---------------------------------------------------------------------------
cat("Reading VCF...\n")
vcf <- tryCatch(
    readVcf(vcf_file, genome = "hg38"),
    error = function(e) {
        cat("ERROR reading VCF:", conditionMessage(e), "\n")
        quit(status = 1)
    }
)
cat("VCF loaded:", nrow(vcf), "variants,", ncol(vcf), "samples\n")

if (nrow(vcf) == 0) {
    cat("WARNING: VCF has 0 variants. Skipping.\n")
    quit(status = 0)
}

# ---------------------------------------------------------------------------
# STEP B: Extract the variant-level table
#
# rowRanges gives us coordinates. We build a clean data frame with:
#   CHROM, POS, ID, REF, ALT, QUAL, FILTER
#
# ALT can be multiple alleles per site (multi-allelic). We collapse them
# to a comma-separated string for easy storage in CSV.
# ---------------------------------------------------------------------------
cat("Extracting variant coordinates...\n")
rr <- rowRanges(vcf)

variant_df <- data.frame(
    CHROM  = as.character(seqnames(rr)),
    POS    = start(rr),
    ID     = as.character(rr$paramRangeID),   # variant ID from VCF
    REF    = as.character(rr$REF),
    ALT    = sapply(rr$ALT, function(a) paste(as.character(a), collapse = ",")),
    QUAL   = rr$QUAL,
    FILTER = sapply(rr$FILTER, function(f) paste(f, collapse = ";")),
    stringsAsFactors = FALSE
)

# Also pull key INFO fields that exist in your VCF
info_data <- info(vcf)
for (field in c("DP", "AC", "AN", "MQ")) {
    if (field %in% names(info_data)) {
        val <- info_data[[field]]
        # AC is a list (one value per ALT allele), collapse it
        if (is.list(val)) {
            variant_df[[field]] <- sapply(val, function(x) paste(x, collapse = ","))
        } else {
            variant_df[[field]] <- as.character(val)
        }
    }
}

cat("Variant table:", nrow(variant_df), "rows,", ncol(variant_df), "columns\n")

# ---------------------------------------------------------------------------
# STEP C: Parse the CSQ (VEP consequence) field
#
# CSQ is the most complex field. Here is how it works:
#
# Each variant can have MULTIPLE CSQ entries separated by commas — one per
# transcript that the variant overlaps. Each entry is pipe-separated (|).
# The column names are defined in the VCF header Description for CSQ.
#
# Example CSQ for one variant:
#   A|missense_variant|MODERATE|CYP2C9|ENSG...|...,
#   A|synonymous_variant|LOW|CYP2C9|ENSG...|...
#
# Strategy: take ONE row per variant — the most severe consequence.
# This avoids the row-count mismatch that crashed your original script
# (which tried to cbind variant table with ALL transcript rows).
#
# Severity order (most to least severe, simplified):
#   stop_gained > frameshift > splice_acceptor > splice_donor >
#   missense_variant > splice_region > synonymous_variant > intronic > ...
# ---------------------------------------------------------------------------
cat("Parsing CSQ field...\n")

# Get CSQ column names from the VCF header
info_hdr <- info(header(vcf))
csq_desc <- as.character(info_hdr["CSQ", "Description"])
csq_format_str <- sub(".*Format: ", "", csq_desc)
csq_cols <- strsplit(csq_format_str, "\\|")[[1]]
cat("CSQ has", length(csq_cols), "fields per transcript\n")

# Define consequence severity (higher index = more severe)
# We use this to pick the most severe transcript consequence per variant
severity_order <- c(
    "intergenic_variant", "downstream_gene_variant", "upstream_gene_variant",
    "non_coding_transcript_variant", "intron_variant", "synonymous_variant",
    "splice_region_variant", "missense_variant", "inframe_insertion",
    "inframe_deletion", "splice_donor_variant", "splice_acceptor_variant",
    "stop_lost", "start_lost", "frameshift_variant", "stop_gained",
    "transcript_ablation"
)

# Function to pick the most severe transcript consequence for one variant
pick_worst_csq <- function(csq_string, csq_cols, severity_order) {

    # If CSQ is missing or empty, return a row of empty strings
    if (is.na(csq_string) || csq_string == "" || csq_string == "character(0)") {
        empty <- rep("", length(csq_cols))
        names(empty) <- csq_cols
        return(empty)
    }

    # Split on comma to get one entry per transcript
    transcripts <- strsplit(csq_string, ",")[[1]]

    # Split each transcript entry on pipe to get field values
    split_entries <- lapply(transcripts, function(t) {
        fields <- strsplit(t, "\\|", fixed = FALSE)[[1]]
        # Pad if fewer fields than expected (truncated entries)
        if (length(fields) < length(csq_cols)) {
            fields <- c(fields, rep("", length(csq_cols) - length(fields)))
        }
        # Trim to expected length if somehow longer
        fields <- fields[seq_len(length(csq_cols))]
        names(fields) <- csq_cols
        return(fields)
    })

    # If only one transcript, return it directly
    if (length(split_entries) == 1) return(split_entries[[1]])

    # Find which consequence column exists
    # VEP names it "Consequence"
    cons_col <- "Consequence"
    if (!cons_col %in% csq_cols) {
        # Fallback: just return first transcript
        return(split_entries[[1]])
    }

    # Score each transcript by severity
    scores <- sapply(split_entries, function(e) {
        cons <- e[cons_col]
        # A variant can have compound consequences like "missense_variant&splice_region"
        # Split on & and take the max severity among them
        parts <- strsplit(cons, "&")[[1]]
        max(match(parts, severity_order, nomatch = 0))
    })

    # Return the transcript with the highest severity score
    return(split_entries[[which.max(scores)]])
}

# Apply to all variants
# info(vcf)$CSQ is a CharacterList — convert to plain character vector
csq_raw <- as.character(info(vcf)$CSQ)

csq_list <- lapply(csq_raw, pick_worst_csq,
                   csq_cols = csq_cols,
                   severity_order = severity_order)

csq_df <- as.data.frame(do.call(rbind, csq_list), stringsAsFactors = FALSE)
colnames(csq_df) <- csq_cols

cat("CSQ table:", nrow(csq_df), "rows (one per variant)\n")

# ---------------------------------------------------------------------------
# STEP D: Extract the genotype matrix
#
# geno(vcf)$GT is a matrix: rows = variants, columns = samples
# Sample names are full BAM paths like:
#   data/AFR/ACB/HG01880/aligned/HG01880.bam
# We clean these to just the sample ID: HG01880
#
# Genotype encoding:
#   0/0 = homozygous reference
#   0/1 = heterozygous
#   1/1 = homozygous alternate
#   ./. = missing
# ---------------------------------------------------------------------------
cat("Extracting genotype matrix...\n")
gt_matrix <- geno(vcf)$GT

# Clean sample names: extract just the sample ID from the BAM path
# Pattern: take the folder name that looks like a sample ID (e.g. HG01880, NA12878)
clean_sample_names <- function(names_vec) {
    sapply(names_vec, function(n) {
        # Try to match a sample ID pattern (letters + digits, 5-7 chars)
        # e.g. HG01880, NA12878, GM12878
        m <- str_extract(n, "[A-Z]{2}\\d{4,6}")
        if (!is.na(m)) return(m)
        # Fallback: take the basename without extension
        return(tools::file_path_sans_ext(basename(n)))
    }, USE.NAMES = FALSE)
}

colnames(gt_matrix) <- clean_sample_names(colnames(gt_matrix))

geno_df <- as.data.frame(gt_matrix, stringsAsFactors = FALSE)
cat("Genotype matrix:", nrow(geno_df), "variants x", ncol(geno_df), "samples\n")

# ---------------------------------------------------------------------------
# STEP E: Save outputs
#
# We save THREE separate files — NOT one merged cbind() table.
# Reason: CSQ and genotype data have different numbers of meaningful columns,
# and cbind() was causing crashes when row counts didn't align perfectly.
#
# Output files:
#   {base}_variants.csv  → CHROM, POS, REF, ALT, QUAL, DP, AC, AN, MQ
#   {base}_csq.csv       → one row per variant, worst-consequence transcript
#   {base}_genotypes.csv → one row per variant, one column per sample (GT calls)
#
# All three files have the same number of rows (one per variant),
# linked by row order — you can merge them later with cbind() safely.
# ---------------------------------------------------------------------------
cat("Saving outputs to:", out_dir, "\n")

write.csv(variant_df,
          file = file.path(out_dir, paste0(base, "_variants.csv")),
          row.names = FALSE)

write.csv(csq_df,
          file = file.path(out_dir, paste0(base, "_csq.csv")),
          row.names = FALSE)

write.csv(geno_df,
          file = file.path(out_dir, paste0(base, "_genotypes.csv")),
          row.names = FALSE)

cat("Done:", base, "\n")
cat("  ->", paste0(base, "_variants.csv"), "\n")
cat("  ->", paste0(base, "_csq.csv"), "\n")
cat("  ->", paste0(base, "_genotypes.csv"), "\n")

REOF

    if [ $? -ne 0 ]; then
        echo "ERROR: R script failed for $VCF"
    else
        echo "SUCCESS: $VCF"
    fi

done

echo ""
echo "All VCFs processed."
