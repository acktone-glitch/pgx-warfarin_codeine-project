# =============================================================================
# 07_loftool_cohort_comparison.R
# Analysis 7: LoFtool gene sensitivity
# Analysis 8: Cohort vs reference AF comparison
# (Combined as both are short analyses)
# =============================================================================

source("D:/RITAH/SCHOOL/Beast Mode/analysis/load_data_full.R")
library(dplyr)
library(ggplot2)

ROOT <- "D:/RITAH/SCHOOL/Beast Mode/analysis"
OUT  <- file.path(ROOT, "results/07_08_loftool_cohort")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# ANALYSIS 7: LoFtool
# =============================================================================
cat("\n============================================\n")
cat("ANALYSIS 7: LoFtool Gene Sensitivity\n")
cat("============================================\n")

loftool_summary <- master %>%
    mutate(LoFtool = as.numeric(LoFtool)) %>%
    filter(!is.na(LoFtool)) %>%
    group_by(Gene, Drug) %>%
    summarise(
        LoFtool_score = round(mean(LoFtool, na.rm = TRUE), 4),
        n_variants    = n(),
        .groups = "drop"
    ) %>%
    arrange(LoFtool_score)

cat("LoFtool scores per gene (lower = more intolerant to LoF):\n")
print(loftool_summary)
write.csv(loftool_summary, file.path(OUT, "loftool_by_gene.csv"), row.names = FALSE)

p1 <- ggplot(loftool_summary,
             aes(x = reorder(Gene, LoFtool_score), y = LoFtool_score, fill = Drug)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    coord_flip() +
    scale_fill_manual(values = c(warfarin = "#2171b5", codeine = "#cb181d")) +
    theme_bw() +
    labs(title = "LoFtool Score per Gene",
         subtitle = "Lower score = gene is more intolerant to loss-of-function variants",
         x = "Gene", y = "LoFtool score",
         caption = "Dashed line = 0.5 midpoint")

ggsave(file.path(OUT, "plot_loftool_by_gene.pdf"), p1, width = 8, height = 5)
ggsave(file.path(OUT, "plot_loftool_by_gene.png"), p1, width = 8, height = 5, dpi = 150)

# =============================================================================
# ANALYSIS 8: Cohort vs Reference AF comparison
# =============================================================================
cat("\n============================================\n")
cat("ANALYSIS 8: Cohort AF vs Reference AF\n")
cat("============================================\n")

# Build matched AF table
cohort_vs_ref <- master %>%
    mutate(
        AF_cohort  = as.numeric(AF_cohort),
        ref_1kg = case_when(
            Population == "AFR" ~ as.numeric(AFR_AF),
            Population == "EUR" ~ as.numeric(EUR_AF),
            Population == "EAS" ~ as.numeric(EAS_AF),
            Population == "SAS" ~ as.numeric(SAS_AF)
        ),
        ref_gnomad = case_when(
            Population == "AFR" ~ as.numeric(gnomADg_AFR_AF),
            Population == "EUR" ~ as.numeric(gnomADg_NFE_AF),
            Population == "EAS" ~ as.numeric(gnomADg_EAS_AF),
            Population == "SAS" ~ as.numeric(gnomADg_SAS_AF)
        )
    ) %>%
    filter(!is.na(AF_cohort)) %>%
    select(Gene, Population, Drug, CHROM, POS, HGVSp,
           AF_cohort, ref_1kg, ref_gnomad,
           AN, Existing_variation, low_AN_flag)

write.csv(cohort_vs_ref, file.path(OUT, "cohort_vs_reference_af.csv"), row.names = FALSE)

# Correlation per population
cat("\nCorrelation between cohort AF and 1KG reference AF:\n")
for (pop in c("AFR","EUR","EAS","SAS")) {
    df <- cohort_vs_ref %>% filter(Population == pop, !is.na(ref_1kg))
    if (nrow(df) > 5) {
        r <- cor(df$AF_cohort, df$ref_1kg, use = "complete.obs", method = "pearson")
        cat(sprintf("  %s: r = %.3f  (n = %d)\n", pop, r, nrow(df)))
    }
}

cat("\nCorrelation between cohort AF and gnomAD genomes AF:\n")
for (pop in c("AFR","EUR","EAS","SAS")) {
    df <- cohort_vs_ref %>% filter(Population == pop, !is.na(ref_gnomad))
    if (nrow(df) > 5) {
        r <- cor(df$AF_cohort, df$ref_gnomad, use = "complete.obs", method = "pearson")
        cat(sprintf("  %s: r = %.3f  (n = %d)\n", pop, r, nrow(df)))
    }
}

# Scatter plot: cohort vs gnomAD
p2 <- ggplot(cohort_vs_ref %>% filter(!is.na(ref_gnomad)),
             aes(x = ref_gnomad, y = AF_cohort, colour = Gene)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~Population) +
    theme_bw() +
    labs(title = "Cohort AF vs gnomAD Genomes AF (matched population)",
         x = "gnomAD genomes AF (matched population)",
         y = "Cohort AF",
         caption = "Dashed line = perfect agreement. Points above = higher in cohort.")

ggsave(file.path(OUT, "plot_cohort_vs_gnomad_scatter.pdf"), p2, width = 10, height = 8)
ggsave(file.path(OUT, "plot_cohort_vs_gnomad_scatter.png"), p2, width = 10, height = 8, dpi = 150)

cat("\nAnalyses 7 & 8 complete. Results saved to:", OUT, "\n")
