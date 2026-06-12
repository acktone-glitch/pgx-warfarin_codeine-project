# =============================================================================
# 01_allele_frequency.R
# Analysis 1: Population allele frequency comparison
#
# QUESTION: Which variants differ in frequency across AFR, EUR, EAS, SAS?
# =============================================================================

source("D:/RITAH/SCHOOL/Beast Mode/analysis/load_data_full.R")
library(dplyr)
library(tidyr)
library(ggplot2)

ROOT <- "D:/RITAH/SCHOOL/Beast Mode/analysis"
OUT  <- file.path(ROOT, "results/01_allele_frequency")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1A. Cohort AF per gene per population
# =============================================================================
cat("\n--- 1A. Cohort AF summary per gene per population ---\n")

af_summary <- master %>%
    mutate(AF_cohort = as.numeric(AF_cohort),
           AN_numeric = as.numeric(AN_numeric)) %>%
    group_by(Gene, Population, Drug) %>%
    summarise(
        n_variants   = n(),
        mean_AF      = round(mean(AF_cohort, na.rm = TRUE), 4),
        median_AF    = round(median(AF_cohort, na.rm = TRUE), 4),
        n_common     = sum(AF_cohort > 0.05, na.rm = TRUE),
        n_rare       = sum(AF_cohort <= 0.05, na.rm = TRUE),
        mean_AN      = round(mean(AN_numeric, na.rm = TRUE), 1),
        .groups = "drop"
    )

print(af_summary, n = Inf)
write.csv(af_summary, file.path(OUT, "af_summary_by_gene_population.csv"), row.names = FALSE)

# =============================================================================
# 1B. Cohort AF vs 1000 Genomes reference AF
# =============================================================================
cat("\n--- 1B. Cohort AF vs 1000 Genomes reference AF ---\n")

af_vs_ref <- master %>%
    mutate(
        AF_cohort  = as.numeric(AF_cohort),
        ref_AF = case_when(
            Population == "AFR" ~ as.numeric(AFR_AF),
            Population == "EUR" ~ as.numeric(EUR_AF),
            Population == "EAS" ~ as.numeric(EAS_AF),
            Population == "SAS" ~ as.numeric(SAS_AF),
            TRUE ~ NA_real_
        )
    ) %>%
    filter(!is.na(ref_AF), !is.na(AF_cohort)) %>%
    select(Gene, Population, CHROM, POS, REF, ALT,
           HGVSp, Consequence, AF_cohort, ref_AF,
           Existing_variation, ClinVar_CLNSIG, low_AN_flag)

cat("Variants with both cohort and 1KG AF:", nrow(af_vs_ref), "\n")
write.csv(af_vs_ref, file.path(OUT, "cohort_vs_1kg_af.csv"), row.names = FALSE)

# =============================================================================
# 1C. Cohort AF vs gnomAD genomes reference AF
# =============================================================================
cat("\n--- 1C. Cohort AF vs gnomAD genomes ---\n")

af_vs_gnomad <- master %>%
    mutate(
        AF_cohort   = as.numeric(AF_cohort),
        gnomad_AF = case_when(
            Population == "AFR" ~ as.numeric(gnomADg_AFR_AF),
            Population == "EUR" ~ as.numeric(gnomADg_NFE_AF),
            Population == "EAS" ~ as.numeric(gnomADg_EAS_AF),
            Population == "SAS" ~ as.numeric(gnomADg_SAS_AF),
            TRUE ~ NA_real_
        )
    ) %>%
    filter(!is.na(gnomad_AF), !is.na(AF_cohort)) %>%
    select(Gene, Population, CHROM, POS, REF, ALT,
           HGVSp, Consequence, AF_cohort, gnomad_AF,
           Existing_variation, ClinVar_CLNSIG, low_AN_flag)

cat("Variants with both cohort and gnomAD AF:", nrow(af_vs_gnomad), "\n")
write.csv(af_vs_gnomad, file.path(OUT, "cohort_vs_gnomad_af.csv"), row.names = FALSE)

# =============================================================================
# 1D. Population-specific variants
#     High in one population (>5%), low in all others (<1%)
# =============================================================================
cat("\n--- 1D. Population-specific variants ---\n")

pop_specific <- master %>%
    mutate(across(c(AFR_AF, EAS_AF, EUR_AF, SAS_AF, AF_cohort), as.numeric)) %>%
    filter(!is.na(AFR_AF), !is.na(EAS_AF), !is.na(EUR_AF), !is.na(SAS_AF)) %>%
    mutate(
        afr_specific = AFR_AF > 0.05 & EAS_AF < 0.01 & EUR_AF < 0.01 & SAS_AF < 0.01,
        eas_specific = EAS_AF > 0.05 & AFR_AF < 0.01 & EUR_AF < 0.01 & SAS_AF < 0.01,
        eur_specific = EUR_AF > 0.05 & AFR_AF < 0.01 & EAS_AF < 0.01 & SAS_AF < 0.01,
        sas_specific = SAS_AF > 0.05 & AFR_AF < 0.01 & EAS_AF < 0.01 & EUR_AF < 0.01
    ) %>%
    filter(afr_specific | eas_specific | eur_specific | sas_specific) %>%
    mutate(specific_to = case_when(
        afr_specific ~ "AFR",
        eas_specific ~ "EAS",
        eur_specific ~ "EUR",
        sas_specific ~ "SAS"
    )) %>%
    select(Gene, specific_to, CHROM, POS, REF, ALT, HGVSp,
           Consequence, IMPACT, AFR_AF, EAS_AF, EUR_AF, SAS_AF,
           AF_cohort, Population, ClinVar_CLNSIG, Existing_variation)

cat("Population-specific variants found:", nrow(pop_specific), "\n")
print(pop_specific %>% select(Gene, specific_to, POS, HGVSp, Consequence, AFR_AF, EAS_AF, EUR_AF, SAS_AF))
write.csv(pop_specific, file.path(OUT, "population_specific_variants.csv"), row.names = FALSE)

# =============================================================================
# 1E. Plots
# =============================================================================
cat("\n--- 1E. Generating plots ---\n")

# Plot 1: Mean AF per gene per population
p1 <- ggplot(af_summary, aes(x = Gene, y = mean_AF, fill = Population)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~Drug, scales = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Mean Allele Frequency per Gene per Population",
         x = "Gene", y = "Mean AF (cohort)", fill = "Population") +
    scale_fill_manual(values = c(AFR="#E41A1C", EUR="#377EB8", EAS="#4DAF4A", SAS="#984EA3"))

ggsave(file.path(OUT, "plot_mean_AF_by_gene_population.pdf"), p1, width = 10, height = 6)
ggsave(file.path(OUT, "plot_mean_AF_by_gene_population.png"), p1, width = 10, height = 6, dpi = 150)

# Plot 2: Cohort AF vs 1KG AF scatter per population
p2 <- ggplot(af_vs_ref, aes(x = ref_AF, y = AF_cohort, colour = Gene)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~Population) +
    theme_bw() +
    labs(title = "Cohort AF vs 1000 Genomes Reference AF",
         x = "1000 Genomes AF (matched population)",
         y = "Cohort AF",
         caption = "Dashed line = perfect agreement") +
    xlim(0, 1) + ylim(0, 1)

ggsave(file.path(OUT, "plot_cohort_vs_1kg_scatter.pdf"), p2, width = 10, height = 8)
ggsave(file.path(OUT, "plot_cohort_vs_1kg_scatter.png"), p2, width = 10, height = 8, dpi = 150)

cat("\nAnalysis 1 complete. Results saved to:", OUT, "\n")
