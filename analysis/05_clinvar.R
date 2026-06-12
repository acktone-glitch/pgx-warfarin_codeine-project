# =============================================================================
# 05_clinvar.R
# Analysis 5: ClinVar drug response variants
#
# QUESTION: Which variants have established clinical evidence for
#           drug response, and how do they differ across populations?
# =============================================================================

source("D:/RITAH/SCHOOL/Beast Mode/analysis/load_data_full.R")
library(dplyr)
library(ggplot2)
library(tidyr)

ROOT <- "D:/RITAH/SCHOOL/Beast Mode/analysis"
OUT  <- file.path(ROOT, "results/05_clinvar")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 5A. All ClinVar variants
# =============================================================================
cat("\n--- 5A. All ClinVar variants ---\n")

clinvar_all <- master %>%
    filter(!is.na(ClinVar_CLNSIG)) %>%
    dplyr::select(Gene, Population, Drug, CHROM, POS, REF, ALT,
           HGVSc, HGVSp, Consequence, IMPACT,
           ClinVar_CLNSIG, ClinVar_CLNDN, ClinVar_CLNREVSTAT,
           AF_cohort, AN, AFR_AF, EAS_AF, EUR_AF, SAS_AF,
           gnomADg_AFR_AF, gnomADg_EAS_AF, gnomADg_NFE_AF, gnomADg_SAS_AF,
           Existing_variation, PUBMED, low_AN_flag)

cat("Variants with ClinVar entries:", nrow(clinvar_all), "\n")
cat("\nClinVar significance breakdown:\n")
print(table(clinvar_all$ClinVar_CLNSIG))

write.csv(clinvar_all, file.path(OUT, "clinvar_all_variants.csv"), row.names = FALSE)

# =============================================================================
# 5B. Drug response variants specifically
# =============================================================================
cat("\n--- 5B. Drug response variants ---\n")

drug_response <- master %>%
    filter(!is.na(ClinVar_CLNSIG) &
           grepl("drug_response", ClinVar_CLNSIG)) %>%
    mutate(AF_cohort = as.numeric(AF_cohort)) %>%
    dplyr::select(Gene, Population, Drug, CHROM, POS, REF, ALT,
           HGVSc, HGVSp, Consequence, IMPACT,
           ClinVar_CLNSIG, ClinVar_CLNDN, ClinVar_CLNREVSTAT,
           AF_cohort, AN, AFR_AF, EAS_AF, EUR_AF, SAS_AF,
           gnomADg_AFR_AF, gnomADg_EAS_AF, gnomADg_NFE_AF, gnomADg_SAS_AF,
           Existing_variation, PUBMED, low_AN_flag) %>%
    arrange(Drug, Gene, Population)

cat("Drug response variants:", nrow(drug_response), "\n")
cat("\nDrug response variants per gene per population:\n")
print(table(drug_response$Gene, drug_response$Population))

cat("\nReview status of drug response entries:\n")
print(table(drug_response$ClinVar_CLNREVSTAT))

write.csv(drug_response, file.path(OUT, "drug_response_variants.csv"), row.names = FALSE)

# =============================================================================
# 5C. Drug response AF across populations
# =============================================================================
cat("\n--- 5C. Drug response variant frequencies ---\n")

dr_freq <- drug_response %>%
    group_by(Gene, POS, REF, ALT, HGVSp, ClinVar_CLNDN) %>%
    summarise(
        AFR_cohort_AF = AF_cohort[Population == "AFR"][1],
        EUR_cohort_AF = AF_cohort[Population == "EUR"][1],
        EAS_cohort_AF = AF_cohort[Population == "EAS"][1],
        SAS_cohort_AF = AF_cohort[Population == "SAS"][1],
        ref_AFR_AF    = as.numeric(AFR_AF)[1],
        ref_EAS_AF    = as.numeric(EAS_AF)[1],
        ref_EUR_AF    = as.numeric(EUR_AF)[1],
        ref_SAS_AF    = as.numeric(SAS_AF)[1],
        .groups = "drop"
    )

print(dr_freq, n = Inf)
write.csv(dr_freq, file.path(OUT, "drug_response_freq_by_population.csv"), row.names = FALSE)

# =============================================================================
# 5D. Plots
# =============================================================================
cat("\n--- 5D. Generating plots ---\n")

# Plot 1: ClinVar significance breakdown per gene
p1 <- clinvar_all %>%
    mutate(CLNSIG_simple = case_when(
        grepl("drug_response", ClinVar_CLNSIG) ~ "drug_response",
        grepl("Pathogenic",    ClinVar_CLNSIG) ~ "Pathogenic",
        grepl("Benign",        ClinVar_CLNSIG) ~ "Benign/Likely_benign",
        grepl("Uncertain",     ClinVar_CLNSIG) ~ "Uncertain_significance",
        TRUE ~ ClinVar_CLNSIG
    )) %>%
    dplyr::count(Gene, CLNSIG_simple) %>%
    ggplot(aes(x = Gene, y = n, fill = CLNSIG_simple)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_bw() +
    scale_fill_manual(values = c(
        "drug_response"          = "#E41A1C",
        "Benign/Likely_benign"   = "#4DAF4A",
        "Uncertain_significance" = "#FF7F00",
        "Pathogenic"             = "#984EA3"
    )) +
    labs(title = "ClinVar Clinical Significance per Gene",
         x = "Gene", y = "Number of variants",
         fill = "ClinVar significance")

ggsave(file.path(OUT, "plot_clinvar_by_gene.pdf"), p1, width = 9, height = 6)
ggsave(file.path(OUT, "plot_clinvar_by_gene.png"), p1, width = 9, height = 6, dpi = 150)

# Plot 2: Drug response AF by population per gene
if (nrow(drug_response) > 0) {
    p2 <- drug_response %>%
        filter(!is.na(AF_cohort)) %>%
        ggplot(aes(x = Population, y = AF_cohort, fill = Population)) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
        facet_wrap(~Gene, scales = "free_y") +
        scale_fill_manual(values = c(AFR="#E41A1C",EUR="#377EB8",EAS="#4DAF4A",SAS="#984EA3")) +
        theme_bw() +
        labs(title = "Drug Response Variant AF by Population",
             x = "Population", y = "Allele Frequency (cohort)")

    ggsave(file.path(OUT, "plot_drug_response_af_by_population.pdf"), p2, width = 12, height = 8)
    ggsave(file.path(OUT, "plot_drug_response_af_by_population.png"), p2, width = 12, height = 8, dpi = 150)
}

cat("\nAnalysis 5 complete. Results saved to:", OUT, "\n")
