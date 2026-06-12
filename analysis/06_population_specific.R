# =============================================================================
# 06_population_specific.R
# Analysis 6: Population-specific variant discovery
#
# QUESTION: Which variants are uniquely enriched in one population?
#           Are any common in AFR but missing from reference databases?
# =============================================================================

source("D:/RITAH/SCHOOL/Beast Mode/analysis/load_data_full.R")
library(dplyr)
library(ggplot2)
library(tidyr)

ROOT <- "D:/RITAH/SCHOOL/Beast Mode/analysis"
OUT  <- file.path(ROOT, "results/06_population_specific")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 6A. Variants absent from gnomAD (novel/rare in reference populations)
# =============================================================================
cat("\n--- 6A. Variants absent or rare in gnomAD ---\n")

novel <- master %>%
    mutate(
        gnomADg_AF = as.numeric(gnomADg_AF),
        AF_cohort  = as.numeric(AF_cohort)
    ) %>%
    filter(gnomADg_AF < 0.001 | is.na(gnomADg_AF)) %>%
    filter(AF_cohort > 0.01) %>%
    dplyr::select(Gene, Population, Drug, CHROM, POS, REF, ALT,
           HGVSc, HGVSp, Consequence, IMPACT,
           AF_cohort, AN, gnomADg_AF, gnomADg_AFR_AF,
           AFR_AF, EAS_AF, EUR_AF, SAS_AF,
           MAX_AF, MAX_AF_POPS, ClinVar_CLNSIG,
           Existing_variation, low_AN_flag) %>%
    arrange(Gene, Population, desc(AF_cohort))

cat("Variants common in cohort (AF>1%) but rare/absent in gnomAD:", nrow(novel), "\n")
print(novel %>% dplyr::select(Gene, Population, HGVSp, Consequence, AF_cohort, gnomADg_AF, ClinVar_CLNSIG))
write.csv(novel, file.path(OUT, "rare_in_gnomad_common_in_cohort.csv"), row.names = FALSE)

# =============================================================================
# 6B. AFR-enriched variants (higher in AFR than all other populations)
# =============================================================================
cat("\n--- 6B. AFR-enriched variants ---\n")

afr_enriched <- master %>%
    mutate(across(c(AFR_AF, EAS_AF, EUR_AF, SAS_AF, AF_cohort,
                    gnomADg_AFR_AF), as.numeric)) %>%
    filter(!is.na(AFR_AF),
           AFR_AF > 0.05,
           (is.na(EAS_AF) | EAS_AF < AFR_AF * 0.5),
           (is.na(EUR_AF) | EUR_AF < AFR_AF * 0.5),
           (is.na(SAS_AF) | SAS_AF < AFR_AF * 0.5)) %>%
    dplyr::select(Gene, Population, CHROM, POS, REF, ALT, HGVSc, HGVSp,
           Consequence, IMPACT, AF_cohort, AFR_AF, EAS_AF, EUR_AF, SAS_AF,
           gnomADg_AFR_AF, ClinVar_CLNSIG, Existing_variation,
           MAX_AF_POPS, low_AN_flag) %>%
    arrange(Gene, desc(AFR_AF))

cat("AFR-enriched variants:", nrow(afr_enriched), "\n")
print(afr_enriched %>% dplyr::select(Gene, HGVSp, Consequence, IMPACT,
                               AFR_AF, EAS_AF, EUR_AF, SAS_AF, ClinVar_CLNSIG))
write.csv(afr_enriched, file.path(OUT, "afr_enriched_variants.csv"), row.names = FALSE)

# =============================================================================
# 6C. MAX_AF_POPS summary — which population has highest AF per variant
# =============================================================================
cat("\n--- 6C. MAX_AF_POPS summary ---\n")

max_pop_summary <- master %>%
    dplyr::count(MAX_AF_POPS, Gene) %>%
    arrange(Gene, desc(n))

cat("Which population has highest AF per variant:\n")
print(table(master$MAX_AF_POPS))
write.csv(max_pop_summary, file.path(OUT, "max_af_pop_by_gene.csv"), row.names = FALSE)

# =============================================================================
# 6D. Plots
# =============================================================================
cat("\n--- 6D. Generating plots ---\n")

# Plot 1: MAX_AF_POPS distribution per gene
p1 <- master %>%
    dplyr::count(Gene, MAX_AF_POPS) %>%
    ggplot(aes(x = Gene, y = n, fill = MAX_AF_POPS)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_bw() +
    labs(title = "Population with Highest Allele Frequency per Variant",
         x = "Gene", y = "Number of variants",
         fill = "Population with\nhighest AF",
         caption = "AFR dominates — reflects greater African genetic diversity")

ggsave(file.path(OUT, "plot_max_af_pops_by_gene.pdf"), p1, width = 10, height = 6)
ggsave(file.path(OUT, "plot_max_af_pops_by_gene.png"), p1, width = 10, height = 6, dpi = 150)

# Plot 2: AFR vs EUR AF scatter — population divergence
ref_data <- master %>%
    mutate(AFR_AF = as.numeric(AFR_AF), EUR_AF = as.numeric(EUR_AF)) %>%
    filter(!is.na(AFR_AF), !is.na(EUR_AF))

if (nrow(ref_data) > 0) {
    p2 <- ggplot(ref_data, aes(x = EUR_AF, y = AFR_AF, colour = Gene)) +
        geom_point(alpha = 0.6, size = 1.8) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
        theme_bw() +
        labs(title = "AFR vs EUR Allele Frequency (1000 Genomes)",
             x = "EUR allele frequency",
             y = "AFR allele frequency",
             caption = "Points above diagonal = higher in AFR; below = higher in EUR")

    ggsave(file.path(OUT, "plot_afr_vs_eur_af.pdf"), p2, width = 8, height = 6)
    ggsave(file.path(OUT, "plot_afr_vs_eur_af.png"), p2, width = 8, height = 6, dpi = 150)
}

cat("\nAnalysis 6 complete. Results saved to:", OUT, "\n")
