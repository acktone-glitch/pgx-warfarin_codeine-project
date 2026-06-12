# =============================================================================
# 04_functional_impact.R
# Analysis 4: Functional impact of missense variants
#
# QUESTION: Which missense variants are likely damaging to protein function?
# =============================================================================

source("D:/RITAH/SCHOOL/Beast Mode/analysis/load_data_full.R")
library(dplyr)
library(ggplot2)
library(tidyr)

ROOT <- "D:/RITAH/SCHOOL/Beast Mode/analysis"
OUT  <- file.path(ROOT, "results/04_functional_impact")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 4A. Extract missense and HIGH impact variants
# =============================================================================
cat("\n--- 4A. Missense and HIGH impact variants ---\n")

missense <- master %>%
    filter(IMPACT %in% c("MODERATE", "HIGH")) %>%
    mutate(
        CADD_phred           = as.numeric(CADD_phred),
        SIFT_score           = suppressWarnings(as.numeric(SIFT_score)),
        Polyphen2_HDIV_score = suppressWarnings(as.numeric(Polyphen2_HDIV_score)),
        REVEL_score          = suppressWarnings(as.numeric(REVEL_score)),
        AF_cohort            = as.numeric(AF_cohort),
        # AlphaMissense — only available for CYP2D6 (dbNSFP5.3)
        AlphaMissense_score  = if ("AlphaMissense_score" %in% names(.))
                                   suppressWarnings(as.numeric(AlphaMissense_score))
                               else NA_real_,
        BayesDel_score       = if ("BayesDel_addAF_score" %in% names(.))
                                   suppressWarnings(as.numeric(BayesDel_addAF_score))
                               else NA_real_
    )

cat("MODERATE/HIGH impact variants:", nrow(missense), "\n")

# =============================================================================
# 4B. Multi-tool scoring — classify each variant
# =============================================================================
cat("\n--- 4B. Multi-tool pathogenicity classification ---\n")

missense <- missense %>%
    mutate(
        SIFT_damaging     = !is.na(SIFT_score) & SIFT_score < 0.05,
        PolyPhen_damaging = !is.na(Polyphen2_HDIV_score) & Polyphen2_HDIV_score > 0.909,
        CADD_high         = !is.na(CADD_phred) & CADD_phred >= 20,
        REVEL_pathogenic  = !is.na(REVEL_score) & REVEL_score > 0.5,

        # AlphaMissense: >0.564 = likely pathogenic (for CYP2D6 only)
        AlphaMissense_damaging = !is.na(AlphaMissense_score) & AlphaMissense_score > 0.564,
        # Count how many tools flag it as damaging
        n_tools_damaging  = SIFT_damaging + PolyPhen_damaging + CADD_high +
                            REVEL_pathogenic + AlphaMissense_damaging,

        # Overall classification
        Pathogenicity_class = case_when(
            n_tools_damaging >= 3 ~ "Likely damaging",
            n_tools_damaging == 2 ~ "Possibly damaging",
            n_tools_damaging == 1 ~ "Uncertain",
            n_tools_damaging == 0 & IMPACT == "HIGH" ~ "HIGH impact (no scores)",
            TRUE ~ "Likely benign / no data"
        )
    )

cat("\nPathogenicity classification:\n")
print(table(missense$Pathogenicity_class))

# =============================================================================
# 4C. Priority variants table
# =============================================================================
priority <- missense %>%
    filter(Pathogenicity_class %in% c("Likely damaging", "Possibly damaging",
                                       "HIGH impact (no scores)")) %>%
    dplyr::select(Gene, Population, Drug, CHROM, POS, REF, ALT,
           HGVSc, HGVSp, Consequence, IMPACT,
           CADD_phred, SIFT_score, Polyphen2_HDIV_score, REVEL_score,
           n_tools_damaging, Pathogenicity_class,
           AF_cohort, AN, AFR_AF, EAS_AF, EUR_AF, SAS_AF,
           DOMAINS, ClinVar_CLNSIG, Existing_variation, low_AN_flag) %>%
    arrange(desc(n_tools_damaging), desc(CADD_phred))

cat("\nPriority variants (damaging by 2+ tools):\n")
print(as.data.frame(priority %>% dplyr::select(Gene, Population, HGVSp, CADD_phred, Pathogenicity_class, ClinVar_CLNSIG, AF_cohort)))

write.csv(missense,  file.path(OUT, "all_missense_scored.csv"),    row.names = FALSE)
write.csv(priority,  file.path(OUT, "priority_damaging_variants.csv"), row.names = FALSE)

# =============================================================================
# 4D. Domain analysis — are damaging variants in functional domains?
# =============================================================================
cat("\n--- 4D. Domain analysis ---\n")

domain_variants <- missense %>%
    filter(!is.na(DOMAINS)) %>%
    dplyr::select(Gene, Population, HGVSp, DOMAINS, CADD_phred,
           Pathogenicity_class, AF_cohort, ClinVar_CLNSIG)

cat("Missense variants with domain annotation:", nrow(domain_variants), "\n")
print(domain_variants)
write.csv(domain_variants, file.path(OUT, "domain_variants.csv"), row.names = FALSE)

# =============================================================================
# 4E. Plots
# =============================================================================
cat("\n--- 4E. Generating plots ---\n")

# CADD score distribution
if (any(!is.na(missense$CADD_phred))) {
    p1 <- ggplot(missense %>% filter(!is.na(CADD_phred)),
                 aes(x = CADD_phred, fill = Gene)) +
        geom_histogram(bins = 20, colour = "white") +
        geom_vline(xintercept = 20, linetype = "dashed", colour = "red") +
        annotate("text", x = 21, y = Inf, label = "CADD=20\n(top 1%)",
                 hjust = 0, vjust = 1.5, colour = "red", size = 3) +
        facet_wrap(~Drug) +
        theme_bw() +
        labs(title = "CADD Score Distribution of Missense Variants",
             x = "CADD phred score", y = "Count",
             caption = "Red line = CADD 20 threshold (top 1% most deleterious)")

    ggsave(file.path(OUT, "plot_cadd_distribution.pdf"), p1, width = 10, height = 5)
    ggsave(file.path(OUT, "plot_cadd_distribution.png"), p1, width = 10, height = 5, dpi = 150)
}

# Pathogenicity class summary
p2 <- ggplot(missense, aes(x = Gene, fill = Pathogenicity_class)) +
    geom_bar() +
    coord_flip() +
    theme_bw() +
    scale_fill_manual(values = c(
        "Likely damaging"          = "#d73027",
        "Possibly damaging"        = "#fc8d59",
        "Uncertain"                = "#fee090",
        "Likely benign / no data"  = "#91bfdb",
        "HIGH impact (no scores)"  = "#7B0000"
    )) +
    labs(title = "Missense Variant Pathogenicity Classification",
         x = "Gene", y = "Number of variants",
         fill = "Classification")

ggsave(file.path(OUT, "plot_pathogenicity_class.pdf"), p2, width = 9, height = 5)
ggsave(file.path(OUT, "plot_pathogenicity_class.png"), p2, width = 9, height = 5, dpi = 150)

cat("\nAnalysis 4 complete. Results saved to:", OUT, "\n")
