# =============================================================================
# 03_star_alleles.R
# Analysis 3: PGx star allele identification
#
# QUESTION: Which known pharmacogenomic star alleles are present,
#           and what are their frequencies across populations?
# =============================================================================

source("D:/RITAH/SCHOOL/Beast Mode/analysis/load_data_full.R")
library(dplyr)
library(ggplot2)
library(tidyr)

ROOT <- "D:/RITAH/SCHOOL/Beast Mode/analysis"
OUT  <- file.path(ROOT, "results/03_star_alleles")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Known PGx star alleles — rsID, HGVSp, gene, allele name, clinical effect
# Source: PharmVar, PharmGKB, CPIC guidelines
# =============================================================================
STAR_ALLELES <- data.frame(
    Gene       = c("CYP2C9","CYP2C9","CYP4F2","VKORC1",
                   "CYP2D6","CYP2D6","CYP2D6","CYP2D6","CYP2D6",
                   "OPRM1","UGT2B7"),
    Star_allele= c("*2","*3","*3","H1(-1639A)",
                   "*4","*5","*6","*10","*17",
                   "118A>G","*2"),
    rsID       = c("rs1799853","rs1057910","rs2108622","rs9923231",
                   "rs3892097","rs5030655","rs5030656","rs1065852","rs28371706",
                   "rs1799971","rs7439366"),
    HGVSp_match= c("p.Arg144Cys","p.Ile359Leu","p.Val433Met",NA,
                   "p.Ile297Leu",NA,NA,"p.Pro34Ser","p.Thr107Ile",
                   "p.Asn40Asp","p.His268Tyr"),
    HGVSc_match= c(NA,NA,NA,"c.-1639G>A",
                   NA,NA,"c.1707delT","c.100C>T","c.320C>T",
                   NA,NA),
    Drug       = c("warfarin","warfarin","warfarin","warfarin",
                   "codeine","codeine","codeine","codeine","codeine",
                   "codeine","codeine"),
    Clinical_effect = c(
        "Reduced CYP2C9 activity — increased warfarin exposure",
        "Severely reduced CYP2C9 activity — high bleeding risk",
        "Reduced vitamin K metabolism — higher warfarin dose needed",
        "Reduced VKORC1 expression — lower warfarin dose needed",
        "Non-functional CYP2D6 — poor metaboliser, no codeine analgesia",
        "Non-functional CYP2D6 — poor metaboliser (gene deletion)",
        "Non-functional CYP2D6 — poor metaboliser (frameshift)",
        "Reduced CYP2D6 activity — common in EAS populations",
        "Reduced CYP2D6 activity — common in AFR populations",
        "Reduced mu-opioid receptor — altered opioid response",
        "Altered morphine glucuronidation ratio"
    ),
    stringsAsFactors = FALSE
)

# CYP2D6 variants in your data — check by position as rsIDs may not match
# due to the different VEP annotation version used
CYP2D6_POSITIONS <- data.frame(
    POS         = c("42129083", "42130522", "42130715"),
    Population  = c("EUR",      "EUR",      "AFR"),
    Notes       = c("EUR variant 1 — check rsID in output",
                    "EUR variant 2 — check rsID in output",
                    "AFR variant (low_MQ_flag=TRUE, interpret cautiously)"),
    stringsAsFactors = FALSE
)
cat("CYP2D6 variants in your data:
")
cyp2d6_in_data <- master[master$Gene == "CYP2D6",
    c("Population","CHROM","POS","REF","ALT","QUAL","MQ",
      "Consequence","HGVSc","HGVSp","Existing_variation",
      "AF_cohort","AN","low_MQ_flag")]
print(cyp2d6_in_data)
cat("\n")

cat("Searching for", nrow(STAR_ALLELES), "known star alleles...\n\n")

# =============================================================================
# Match by rsID (Existing_variation), HGVSp, or HGVSc
# =============================================================================
results <- list()

for (i in seq_len(nrow(STAR_ALLELES))) {
    sa <- STAR_ALLELES[i, ]

    # Match by rsID
    matches <- master[!is.na(master$Existing_variation) &
                      grepl(sa$rsID, master$Existing_variation, fixed = TRUE), ]

    # Also try HGVSp match if rsID not found
    if (nrow(matches) == 0 && !is.na(sa$HGVSp_match)) {
        matches <- master[!is.na(master$HGVSp) &
                          grepl(sa$HGVSp_match, master$HGVSp, fixed = TRUE), ]
    }

    # Also try HGVSc match
    if (nrow(matches) == 0 && !is.na(sa$HGVSc_match)) {
        matches <- master[!is.na(master$HGVSc) &
                          grepl(sa$HGVSc_match, master$HGVSc, fixed = TRUE), ]
    }

    if (nrow(matches) > 0) {
        matches$Star_allele     <- sa$Star_allele
        matches$Clinical_effect <- sa$Clinical_effect
        results[[i]] <- matches %>%
            dplyr::select(Gene, Star_allele, Population, CHROM, POS, REF, ALT,
                   HGVSc, HGVSp, Existing_variation, AF_cohort, AN,
                   AFR_AF, EAS_AF, EUR_AF, SAS_AF,
                   gnomADg_AFR_AF, gnomADg_EAS_AF, gnomADg_NFE_AF, gnomADg_SAS_AF,
                   ClinVar_CLNSIG, Clinical_effect, low_AN_flag)
        cat("FOUND:", sa$Gene, sa$Star_allele, "(", sa$rsID, ") —",
            nrow(matches), "population row(s)\n")
    } else {
        cat("NOT FOUND:", sa$Gene, sa$Star_allele, "(", sa$rsID, ")\n")
    }
}

# =============================================================================
# Combine and summarise
# =============================================================================
if (length(results) > 0) {
    star_found <- bind_rows(results[!sapply(results, is.null)])

    cat("\n============================================\n")
    cat("STAR ALLELES IDENTIFIED:", nrow(star_found), "rows\n")
    cat("============================================\n")
    print_df <- as.data.frame(star_found %>%
        dplyr::select(Gene, Star_allele, Population, AF_cohort, AN,
               HGVSp, ClinVar_CLNSIG, Clinical_effect) %>%
        arrange(Gene, Star_allele, Population))
    print(print_df)

    # Frequency table: star allele x population
    freq_table <- star_found %>%
        mutate(AF_cohort = as.numeric(AF_cohort)) %>%
        dplyr::select(Gene, Star_allele, Population, AF_cohort, AN, Clinical_effect) %>%
        pivot_wider(names_from = Population,
                    values_from = c(AF_cohort, AN),
                    values_fill = list(AF_cohort = NA, AN = 0))

    write.csv(star_found,  file.path(OUT, "star_alleles_found.csv"),   row.names = FALSE)
    write.csv(freq_table,  file.path(OUT, "star_allele_frequencies.csv"), row.names = FALSE)

    # Plot: star allele frequencies across populations
    plot_data <- star_found %>%
        mutate(AF_cohort = as.numeric(AF_cohort)) %>%
        filter(!is.na(AF_cohort))

    if (nrow(plot_data) > 0) {
        plot_data$AF_cohort <- as.numeric(plot_data$AF_cohort)
        p1 <- ggplot(plot_data,
                     aes(x = Population, y = AF_cohort,
                         fill = Population, label = round(AF_cohort, 3))) +
            geom_bar(stat = "identity") +
            geom_text(vjust = -0.3, size = 3) +
            facet_wrap(~paste(Gene, Star_allele), scales = "free_y") +
            scale_fill_manual(values = c(AFR="#E41A1C",EUR="#377EB8",EAS="#4DAF4A",SAS="#984EA3")) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(title = "PGx Star Allele Frequencies by Population",
                 x = "Population", y = "Allele Frequency (cohort)",
                 caption = "Frequency based on cohort AC/AN")

        ggsave(file.path(OUT, "plot_star_allele_frequencies.pdf"), p1, width = 12, height = 8)
        ggsave(file.path(OUT, "plot_star_allele_frequencies.png"), p1, width = 12, height = 8, dpi = 150)
    }
} else {
    cat("\nNo star alleles identified — check rsID matches manually.\n")
}

cat("\nAnalysis 3 complete. Results saved to:", OUT, "\n")
