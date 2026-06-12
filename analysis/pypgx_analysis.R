# =============================================================================
#  pypgx_analysis.R
#  Analyse PyPGx diplotype / phenotype / genotype results
#
#  INPUT:  TSV files from pypgx print-data (one per gene per population)
#          Located at: BASE_DIR/<POP>/<GENE>/<GENE>_<POP>_results.tsv
#
#  OUTPUT: results/pypgx_analysis/ — plots + tables
#
#  USAGE:  Source in RStudio after setting BASE_DIR
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(stringr)

# ── PATHS ─────────────────────────────────────────────────────────────────────
BASE_DIR <- "D:/RITAH/SCHOOL/Beast Mode/analysis/analysis2/pypgx"
OUT_DIR  <- "D:/RITAH/SCHOOL/Beast Mode/analysis/analysis2/results/pypgx_analysis"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

POPS  <- c("AFR", "EAS", "EUR", "SAS")
GENES <- c("CYP2C9", "CYP2C19", "CYP2D6", "CYP4F2", "VKORC1", "UGT2B7", "OPRM1")

POP_COLORS  <- c(AFR="#D62728", EAS="#2CA02C", EUR="#1F77B4", SAS="#9467BD")
GENE_COLORS <- c(CYP2C9="#377EB8", CYP2C19="#984EA3", CYP2D6="#4DAF4A",
                 CYP4F2="#FF7F00", VKORC1="#00CED1", UGT2B7="#999999",
                 OPRM1="#F781BF")

DRUG_MAP <- c(CYP2C9="warfarin", CYP2C19="warfarin", CYP4F2="warfarin",
              VKORC1="warfarin", CYP2D6="codeine", UGT2B7="codeine",
              OPRM1="codeine")

theme_pgx <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(strip.background  = element_rect(fill = "grey92", colour = "grey70"),
          strip.text        = element_text(face = "bold"),
          panel.grid.major  = element_line(colour = "grey90"),
          panel.grid.minor  = element_blank(),
          plot.title        = element_text(face = "bold", size = base_size + 1),
          plot.caption      = element_text(colour = "grey50", size = base_size - 3),
          axis.text         = element_text(colour = "grey20"))
}

save_plot <- function(p, name, w = 10, h = 6) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  message("  Saved: ", name)
}

# =============================================================================
#  1. LOAD ALL RESULTS
# =============================================================================
message("=== Loading PyPGx results ===")

all_results <- list()

for (pop in POPS) {
  for (gene in GENES) {
    tsv <- file.path(BASE_DIR, pop, gene,
                     paste0(gene, "_", pop, "_results.tsv"))
    if (!file.exists(tsv)) {
      message("  MISSING: ", pop, "/", gene)
      next
    }
    df <- tryCatch(
      read_tsv(tsv, show_col_types = FALSE),
      error = function(e) { message("  ERROR reading ", tsv, ": ", e$message); NULL }
    )
    if (is.null(df) || nrow(df) == 0) next

    df$Population <- pop
    df$Gene       <- gene
    df$Drug       <- DRUG_MAP[gene]
    all_results[[paste(pop, gene, sep = "_")]] <- df
    message("  Loaded: ", pop, "/", gene, " (", nrow(df), " samples)")
  }
}

if (length(all_results) == 0) stop("No results loaded. Check BASE_DIR and TSV paths.")

master_pgx <- bind_rows(all_results)
message("\nTotal sample-gene rows: ", nrow(master_pgx))
message("Columns: ", paste(names(master_pgx), collapse = ", "))

# Save combined table
write_csv(master_pgx, file.path(OUT_DIR, "pypgx_all_results.csv"))

# ── Identify column names (vary between genes) ───────────────────────────────
dip_col   <- names(master_pgx)[str_detect(names(master_pgx),
               regex("diplotype|genotype", ignore_case = TRUE))][1]
pheno_col <- names(master_pgx)[str_detect(names(master_pgx),
               regex("phenotype", ignore_case = TRUE))][1]
score_col <- names(master_pgx)[str_detect(names(master_pgx),
               regex("activity.*score|score", ignore_case = TRUE))][1]

message("Diplotype column:  ", dip_col)
message("Phenotype column:  ", pheno_col)
message("Activity score:    ", score_col)

# =============================================================================
#  2. DIPLOTYPE FREQUENCY TABLE
# =============================================================================
message("\n=== Analysis: Diplotype frequencies ===")

if (!is.na(dip_col)) {
  dip_freq <- master_pgx %>%
    rename(Diplotype = all_of(dip_col)) %>%
    filter(!is.na(Diplotype), Diplotype != "Unresolved") %>%
    count(Gene, Population, Diplotype) %>%
    group_by(Gene, Population) %>%
    mutate(Freq = n / sum(n)) %>%
    ungroup() %>%
    arrange(Gene, Population, desc(Freq))

  write_csv(dip_freq, file.path(OUT_DIR, "diplotype_frequencies.csv"))
  message("  Diplotype table: ", nrow(dip_freq), " rows")

  # Top 5 diplotypes per gene per population
  top_dips <- dip_freq %>%
    group_by(Gene, Population) %>%
    slice_max(Freq, n = 5) %>%
    ungroup()
  write_csv(top_dips, file.path(OUT_DIR, "top5_diplotypes_per_gene_pop.csv"))
}

# =============================================================================
#  3. PHENOTYPE DISTRIBUTION
# =============================================================================
message("\n=== Analysis: Phenotype distribution ===")

if (!is.na(pheno_col)) {
  pheno_df <- master_pgx %>%
    rename(Phenotype = all_of(pheno_col)) %>%
    filter(!is.na(Phenotype))

  pheno_count <- pheno_df %>%
    count(Gene, Population, Drug, Phenotype) %>%
    group_by(Gene, Population) %>%
    mutate(Pct = 100 * n / sum(n)) %>%
    ungroup()

  write_csv(pheno_count, file.path(OUT_DIR, "phenotype_counts.csv"))

  # ── Plot 1: Phenotype % stacked bar per gene per population ────────────────
  # Define phenotype colour scale (CPIC standard colours)
  pheno_levels <- c("Poor Metabolizer", "Intermediate Metabolizer",
                    "Normal Metabolizer", "Rapid Metabolizer",
                    "Ultrarapid Metabolizer", "Indeterminate",
                    "Decreased Sensitivity", "Increased Sensitivity",
                    "Normal Sensitivity")

  pheno_colors <- c(
    "Poor Metabolizer"         = "#B2182B",
    "Intermediate Metabolizer" = "#EF8A62",
    "Normal Metabolizer"       = "#4DAF4A",
    "Rapid Metabolizer"        = "#2CA02C",
    "Ultrarapid Metabolizer"   = "#1F77B4",
    "Decreased Sensitivity"    = "#FF7F00",
    "Normal Sensitivity"       = "#984EA3",
    "Increased Sensitivity"    = "#00CED1",
    "Indeterminate"            = "#AAAAAA"
  )

  p_pheno_stacked <- pheno_count %>%
    mutate(Phenotype = factor(Phenotype,
             levels = c(pheno_levels,
                        setdiff(unique(Phenotype), pheno_levels)))) %>%
    ggplot(aes(x = Population, y = Pct, fill = Phenotype)) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = pheno_colors, na.value = "grey70",
                      name = "Phenotype") +
    theme_pgx() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.text = element_text(size = 8)) +
    labs(title = "PyPGx Predicted Phenotype Distribution by Gene and Population",
         x = "Population", y = "Percentage of samples (%)",
         caption = "Based on PyPGx diplotype-to-phenotype translation (CPIC guidelines)")
  save_plot(p_pheno_stacked, "plot_phenotype_distribution_stacked", w = 14, h = 10)

  # ── Plot 2: Poor/Intermediate Metabolizer % per gene per population ─────────
  pm_im <- pheno_count %>%
    filter(str_detect(Phenotype, "Poor|Intermediate|Decreased")) %>%
    group_by(Gene, Population, Drug) %>%
    summarise(Pct_reduced = sum(Pct), .groups = "drop")

  if (nrow(pm_im) > 0) {
    p_pm_im <- ggplot(pm_im, aes(x = Population, y = Pct_reduced, fill = Population)) +
      geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
      geom_text(aes(label = round(Pct_reduced, 1)),
                vjust = -0.3, size = 2.8, colour = "grey20") +
      facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
      scale_fill_manual(values = POP_COLORS) +
      theme_pgx() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(title = "% Samples with Reduced Metaboliser Function (Poor + Intermediate)",
           x = "Population", y = "Percentage (%)",
           caption = "Poor Metabolizer + Intermediate Metabolizer combined")
    save_plot(p_pm_im, "plot_reduced_metabolizer_pct", w = 14, h = 8)
  }

  # ── Summary table: phenotype % per gene per pop ────────────────────────────
  pheno_wide <- pheno_count %>%
    select(Gene, Population, Phenotype, Pct) %>%
    mutate(Pct = round(Pct, 1)) %>%
    pivot_wider(names_from = Population, values_from = Pct, values_fill = 0)
  write_csv(pheno_wide, file.path(OUT_DIR, "phenotype_pct_wide.csv"))
}

# =============================================================================
#  4. ACTIVITY SCORE DISTRIBUTION
# =============================================================================
message("\n=== Analysis: Activity scores ===")

if (!is.na(score_col)) {
  score_df <- master_pgx %>%
    rename(ActivityScore = all_of(score_col)) %>%
    mutate(ActivityScore = suppressWarnings(as.numeric(ActivityScore))) %>%
    filter(!is.na(ActivityScore))

  if (nrow(score_df) > 0) {
    p_score <- ggplot(score_df,
                      aes(x = ActivityScore, fill = Population)) +
      geom_histogram(bins = 15, colour = "white", linewidth = 0.3,
                     position = "dodge") +
      facet_wrap(~ Gene, scales = "free", ncol = 4) +
      scale_fill_manual(values = POP_COLORS) +
      theme_pgx() +
      labs(title = "Activity Score Distribution by Gene and Population (PyPGx)",
           x = "Activity Score", y = "Count",
           caption = "Activity score: 0 = no function, 1 = normal function per allele")
    save_plot(p_score, "plot_activity_score_distribution", w = 14, h = 8)

    # Median activity score per gene per population
    score_summary <- score_df %>%
      group_by(Gene, Population, Drug) %>%
      summarise(median_score = median(ActivityScore),
                mean_score   = round(mean(ActivityScore), 2),
                n_zero       = sum(ActivityScore == 0),
                n_reduced    = sum(ActivityScore < 1),
                n_total      = n(),
                .groups = "drop")
    write_csv(score_summary, file.path(OUT_DIR, "activity_score_summary.csv"))
    message("  Activity score summary:")
    print(score_summary)
  }
}

# =============================================================================
#  5. STAR ALLELE FREQUENCY TABLE
# =============================================================================
message("\n=== Analysis: Star allele frequencies ===")

if (!is.na(dip_col)) {
  # Split diplotype into two haplotypes and count individual allele frequencies
  allele_freq <- master_pgx %>%
    rename(Diplotype = all_of(dip_col)) %>%
    filter(!is.na(Diplotype), Diplotype != "Unresolved",
           str_detect(Diplotype, "/")) %>%
    mutate(
      Hap1 = str_extract(Diplotype, "^[^/]+"),
      Hap2 = str_extract(Diplotype, "[^/]+$")
    ) %>%
    pivot_longer(cols = c(Hap1, Hap2), values_to = "Allele") %>%
    filter(!is.na(Allele)) %>%
    count(Gene, Population, Allele) %>%
    group_by(Gene, Population) %>%
    mutate(AF = round(n / sum(n), 3)) %>%
    ungroup() %>%
    arrange(Gene, Population, desc(AF))

  write_csv(allele_freq, file.path(OUT_DIR, "star_allele_frequencies_pypgx.csv"))
  message("  Star allele frequency table: ", nrow(allele_freq), " rows")

  # Plot: top alleles per gene
  top_alleles <- allele_freq %>%
    group_by(Gene) %>%
    filter(Allele != paste0(Gene, "*1")) %>%   # exclude reference allele
    slice_max(AF, n = 20) %>%
    ungroup()

  if (nrow(top_alleles) > 0) {
    p_allele <- ggplot(top_alleles,
                       aes(x = reorder(Allele, AF), y = AF, fill = Population)) +
      geom_bar(stat = "identity", position = "dodge",
               colour = "white", linewidth = 0.3) +
      coord_flip() +
      facet_wrap(~ Gene, scales = "free", ncol = 3) +
      scale_fill_manual(values = POP_COLORS) +
      theme_pgx() +
      theme(axis.text.y = element_text(size = 7)) +
      labs(title = "Non-Reference Star Allele Frequencies by Population (PyPGx)",
           x = "Star Allele", y = "Allele Frequency",
           caption = "Reference (*1) allele excluded")
    save_plot(p_allele, "plot_star_allele_freq_pypgx", w = 16, h = 12)
  }
}

# =============================================================================
#  6. CROSS-POPULATION COMPARISON: PM/IM HEATMAP
# =============================================================================
message("\n=== Analysis: Reduced function heatmap ===")

if (!is.na(pheno_col)) {
  pm_heat <- master_pgx %>%
    rename(Phenotype = all_of(pheno_col)) %>%
    filter(!is.na(Phenotype)) %>%
    group_by(Gene, Population) %>%
    summarise(
      Pct_reduced = 100 * mean(str_detect(Phenotype,
                    "Poor|Intermediate|Decreased"), na.rm = TRUE),
      .groups = "drop"
    )

  p_heat <- ggplot(pm_heat,
                   aes(x = Population, y = Gene, fill = Pct_reduced)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = round(Pct_reduced, 1)),
              size = 3.5, colour = "grey10") +
    scale_fill_gradient(low = "#EBF3FA", high = "#B2182B",
                        name = "% Reduced\nfunction") +
    theme_pgx() +
    theme(axis.text.x = element_text(angle = 0),
          axis.text.y = element_text(face = "bold")) +
    labs(title = "Proportion of Samples with Reduced Metaboliser Function",
         subtitle = "Poor Metabolizer + Intermediate Metabolizer",
         x = "Population", y = "Gene",
         caption = "PyPGx diplotype-to-phenotype translation")
  save_plot(p_heat, "plot_reduced_function_heatmap", w = 8, h = 6)
}

# =============================================================================
#  7. WARFARIN / CODEINE COMBINED METABOLISER PROFILE
# =============================================================================
message("\n=== Analysis: Drug pathway profiles ===")

if (!is.na(pheno_col)) {
  # Warfarin: combine CYP2C9 + VKORC1 + CYP4F2 per sample
  # Here we summarise at population level since we don't have sample IDs linked
  pathway_summary <- master_pgx %>%
    rename(Phenotype = all_of(pheno_col)) %>%
    filter(!is.na(Phenotype), !is.na(Drug)) %>%
    count(Drug, Gene, Population, Phenotype) %>%
    group_by(Drug, Gene, Population) %>%
    mutate(Pct = 100 * n / sum(n)) %>%
    ungroup()

  p_pathway <- pathway_summary %>%
    mutate(Phenotype = factor(Phenotype)) %>%
    filter(str_detect(Phenotype, "Poor|Intermediate|Normal|Rapid|Ultra|Decreased|Increased")) %>%
    ggplot(aes(x = Population, y = Pct, fill = Phenotype)) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
    facet_grid(Drug ~ Gene, scales = "free_x") +
    scale_fill_manual(values = c(
      "Poor Metabolizer"         = "#B2182B",
      "Intermediate Metabolizer" = "#EF8A62",
      "Normal Metabolizer"       = "#4DAF4A",
      "Rapid Metabolizer"        = "#2CA02C",
      "Ultrarapid Metabolizer"   = "#1F77B4",
      "Decreased Sensitivity"    = "#FF7F00",
      "Normal Sensitivity"       = "#984EA3",
      "Increased Sensitivity"    = "#00CED1"
    ), na.value = "grey70") +
    theme_pgx() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.text = element_text(size = 7)) +
    labs(title = "Warfarin and Codeine Metaboliser Phenotype Profiles by Population",
         x = "Population", y = "Percentage (%)",
         caption = "PyPGx diplotype-to-phenotype translation per gene")
  save_plot(p_pathway, "plot_drug_pathway_profiles", w = 16, h = 8)
}

# =============================================================================
#  SUMMARY
# =============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("PyPGx ANALYSIS COMPLETE")
message(paste(rep("=", 60), collapse = ""))
message("Output: ", OUT_DIR)
message("Plots saved:")
message("  plot_phenotype_distribution_stacked")
message("  plot_reduced_metabolizer_pct")
message("  plot_activity_score_distribution")
message("  plot_star_allele_freq_pypgx")
message("  plot_reduced_function_heatmap")
message("  plot_drug_pathway_profiles")
message("Tables saved:")
message("  pypgx_all_results.csv")
message("  diplotype_frequencies.csv")
message("  phenotype_counts.csv")
message("  phenotype_pct_wide.csv")
message("  activity_score_summary.csv")
message("  star_allele_frequencies_pypgx.csv")
