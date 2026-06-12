#!/usr/bin/env Rscript
# =============================================================================
#  parse_vep_gatk.R
#  Parse VEP-annotated GATK joint-called VCFs into master analysis table
#
#  INPUT:  results/vep/<POP>/<POP>_pgx_filtered.vep.vcf.gz  (4 populations)
#  OUTPUT: results/master_gatk.csv  — combined table ready for R analyses
#
#  USAGE (local RStudio):
#    source("parse_vep_gatk.R")
#  USAGE (HPC):
#    Rscript parse_vep_gatk.R
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

BASE  <- "/etc/ace-data/home/rnagaba/reset"
OUT   <- file.path(BASE, "results")
POPS  <- c("AFR", "EAS", "EUR", "SAS")

# =============================================================================
#  HELPER: parse one VEP VCF into a flat data frame
# =============================================================================
parse_vep_vcf <- function(vcf_path, population) {

  message("  Reading: ", basename(vcf_path))

  # ── Read raw lines ──────────────────────────────────────────────────────────
  con   <- gzcon(file(vcf_path, "rb"))
  lines <- readLines(con)
  close(con)

  # ── Extract CSQ field order from VEP header ─────────────────────────────────
  csq_header <- grep("ID=CSQ", lines, value = TRUE)
  if (length(csq_header) == 0) stop("No CSQ header found in ", vcf_path)

  csq_fields <- csq_header %>%
    str_extract("Format: ([^\"]+)") %>%
    str_replace("Format: ", "") %>%
    str_split("\\|") %>%
    .[[1]]

  message("  CSQ fields: ", length(csq_fields))

  # ── Extract VCF column header ───────────────────────────────────────────────
  col_header_line <- grep("^#CHROM", lines, value = TRUE)
  col_names <- str_split(col_header_line, "\t")[[1]]
  col_names[1] <- "CHROM"

  sample_cols <- col_names[10:length(col_names)]
  n_samples   <- length(sample_cols)
  message("  Samples: ", n_samples)

  # ── Extract data lines ──────────────────────────────────────────────────────
  data_lines <- lines[!str_starts(lines, "#")]
  if (length(data_lines) == 0) {
    warning("No variants found in ", vcf_path)
    return(NULL)
  }

  message("  Variants (raw): ", length(data_lines))

  # ── Parse each line ─────────────────────────────────────────────────────────
  rows <- lapply(data_lines, function(line) {

    fields <- str_split(line, "\t")[[1]]
    if (length(fields) < 8) return(NULL)

    CHROM  <- fields[1]
    POS    <- fields[2]
    ID     <- fields[3]
    REF    <- fields[4]
    ALT    <- fields[5]
    QUAL   <- fields[6]
    FILTER <- fields[7]
    INFO   <- fields[8]
    FORMAT <- if (length(fields) >= 9) fields[9] else NA

    # ── Extract AC, AN, AF from INFO ──────────────────────────────────────────
    AC <- str_extract(INFO, "(?<=AC=)[^;]+") %>% str_split(",") %>% .[[1]] %>% .[1]
    AN <- str_extract(INFO, "(?<=AN=)[^;]+")
    AF <- if (!is.na(AN) && as.numeric(AN) > 0)
            as.character(round(as.numeric(AC) / as.numeric(AN), 6))
          else NA

    # ── Extract CSQ annotations ───────────────────────────────────────────────
    csq_raw <- str_extract(INFO, "(?<=CSQ=)[^\\t]+")
    if (is.na(csq_raw)) return(NULL)

    # Take first CSQ transcript (primary annotation)
    csq_entries <- str_split(csq_raw, ",")[[1]]

    lapply(csq_entries, function(entry) {
      vals <- str_split(entry, "\\|")[[1]]
      if (length(vals) < length(csq_fields)) {
        vals <- c(vals, rep(NA, length(csq_fields) - length(vals)))
      }
      csq_df <- setNames(as.list(vals[seq_along(csq_fields)]), csq_fields)

      c(
        list(
          CHROM      = CHROM,
          POS        = POS,
          ID         = ID,
          REF        = REF,
          ALT        = ALT,
          QUAL       = QUAL,
          FILTER     = FILTER,
          AC         = AC,
          AN         = AN,
          AF_cohort  = AF,
          Population = population
        ),
        csq_df
      )
    })
  })

  # ── Flatten nested list ─────────────────────────────────────────────────────
  rows_flat <- Filter(Negate(is.null), unlist(rows, recursive = FALSE))

  if (length(rows_flat) == 0) return(NULL)

  df <- bind_rows(lapply(rows_flat, as.data.frame, stringsAsFactors = FALSE))

  message("  Rows after CSQ expansion: ", nrow(df))
  df
}

# =============================================================================
#  PARSE ALL 4 POPULATIONS
# =============================================================================
all_pops <- list()

for (pop in POPS) {
  vcf_path <- file.path(BASE, "results", "vep", pop,
                        paste0(pop, "_pgx_filtered.vep.vcf.gz"))

  if (!file.exists(vcf_path)) {
    warning("VCF not found — skipping: ", vcf_path)
    next
  }

  message("\n=== ", pop, " ===")
  df <- tryCatch(
    parse_vep_vcf(vcf_path, pop),
    error = function(e) { warning(pop, ": ", e$message); NULL }
  )

  if (!is.null(df)) all_pops[[pop]] <- df
}

if (length(all_pops) == 0) stop("No populations parsed successfully.")

# =============================================================================
#  COMBINE AND CLEAN
# =============================================================================
message("\n=== Combining populations ===")
master <- bind_rows(all_pops)

message("Total rows before filtering: ", nrow(master))

# ── Keep MANE Select or canonical transcript only ──────────────────────────
# Priority: MANE_SELECT > CANONICAL > everything else
# This avoids the alternate transcript problem from the previous analysis
if ("MANE_SELECT" %in% names(master)) {
  master <- master %>%
    group_by(CHROM, POS, REF, ALT, Population) %>%
    arrange(
      desc(!is.na(MANE_SELECT) & MANE_SELECT != ""),
      desc(CANONICAL == "YES")
    ) %>%
    slice(1) %>%
    ungroup()
  message("Rows after MANE Select / canonical filter: ", nrow(master))
}

# ── Numeric conversions ────────────────────────────────────────────────────
numeric_cols <- c("AC", "AN", "AF_cohort", "CADD_phred", "REVEL_score",
                  "AlphaMissense_score", "GERP++_RS", "phyloP100way_vertebrate",
                  "gnomAD4.1_joint_AF", "gnomAD4.1_joint_AFR_AF",
                  "gnomAD4.1_joint_EAS_AF", "gnomAD4.1_joint_NFE_AF",
                  "gnomAD4.1_joint_SAS_AF",
                  "gnomAD2.1.1_exomes_controls_AF",
                  "1000Gp3_AF", "1000Gp3_AFR_AF",
                  "1000Gp3_EAS_AF", "1000Gp3_EUR_AF", "1000Gp3_SAS_AF")

for (col in numeric_cols) {
  if (col %in% names(master)) {
    master[[col]] <- suppressWarnings(as.numeric(master[[col]]))
  }
}

# ── Replace empty strings with NA ─────────────────────────────────────────
master <- master %>%
  mutate(across(where(is.character), ~ na_if(., "")))

# ── Add drug pathway annotation ────────────────────────────────────────────
warfarin_genes <- c("CYP2C9", "CYP2C19", "VKORC1", "CYP4F2", "GGCX", "CALU")
codeine_genes  <- c("CYP2D6", "UGT2B7", "OPRM1", "DRD2")

master <- master %>%
  mutate(Drug = case_when(
    SYMBOL %in% warfarin_genes ~ "warfarin",
    SYMBOL %in% codeine_genes  ~ "codeine",
    TRUE                       ~ "other"
  ))

# ── Add AN quality flag ─────────────────────────────────────────────────────
master <- master %>%
  mutate(
    AN_numeric   = suppressWarnings(as.numeric(AN)),
    low_AN_flag  = !is.na(AN_numeric) & AN_numeric < 4,
    MQ_flag      = if ("MQ" %in% names(.)) {
                     suppressWarnings(as.numeric(MQ)) < 30
                   } else FALSE
  )

# ── FILTER: keep only PASS variants ────────────────────────────────────────
message("PASS variants: ", sum(master$FILTER == "PASS", na.rm = TRUE),
        " of ", nrow(master))
master_pass <- master %>% filter(FILTER == "PASS")
message("Rows after PASS filter: ", nrow(master_pass))

# =============================================================================
#  SAVE OUTPUTS
# =============================================================================
dir.create(file.path(OUT), showWarnings = FALSE)

# Full table (all variants including non-PASS — for QC)
out_all  <- file.path(OUT, "master_gatk_all.csv")
# PASS-only table (for analysis)
out_pass <- file.path(OUT, "master_gatk_pass.csv")

write_csv(master,      out_all)
write_csv(master_pass, out_pass)

message("\n=== DONE ===")
message("All variants:   ", nrow(master),      " rows → ", out_all)
message("PASS variants:  ", nrow(master_pass), " rows → ", out_pass)
message("Populations:    ", paste(unique(master$Population), collapse = ", "))

if ("SYMBOL" %in% names(master_pass)) {
  gene_counts <- table(master_pass$SYMBOL, master_pass$Population)
  message("\nVariants per gene per population (PASS, MANE Select):")
  print(gene_counts)
}
