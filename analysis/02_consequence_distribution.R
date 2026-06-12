# =============================================================================
# 02_consequence_distribution.R
# Analysis 2: Consequence and IMPACT distribution per gene and population
#
# QUESTION: What types of variants are present and where are they concentrated?
# =============================================================================

source("D:/RITAH/SCHOOL/Beast Mode/analysis/load_data_full.R")
library(dplyr)
library(ggplot2)
library(tidyr)

ROOT <- "D:/RITAH/SCHOOL/Beast Mode/analysis"
OUT  <- file.path(ROOT, "results/02_consequence")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2A. Overall consequence counts
# =============================================================================
cat("\n--- 2A. Overall consequence counts ---\n")

# Simplify compound consequences (e.g. "intron_variant&non_coding_transcript_variant")
master$Consequence_simple <- sapply(master$Consequence, function(x) {
    if (is.na(x)) return(NA)
    strsplit(x, "&")[[1]][1]   # take most severe (first) component
})

cons_overall <- master %>%
    dplyr::count(Consequence_simple, IMPACT, sort = TRUE) %>%
    dplyr::rename(Consequence = Consequence_simple, N = n)

print(cons_overall)
write.csv(cons_overall, file.path(OUT, "consequence_overall.csv"), row.names = FALSE)

# =============================================================================
# 2B. IMPACT breakdown per gene
# =============================================================================
cat("\n--- 2B. IMPACT per gene ---\n")

impact_by_gene <- master %>%
    dplyr::count(Gene, Drug, IMPACT) %>%
    pivot_wider(names_from = IMPACT, values_from = n, values_fill = 0) %>%
    arrange(Drug, Gene)

print(impact_by_gene)
write.csv(impact_by_gene, file.path(OUT, "impact_by_gene.csv"), row.names = FALSE)

# =============================================================================
# 2C. IMPACT breakdown per population
# =============================================================================
cat("\n--- 2C. IMPACT per population ---\n")

impact_by_pop <- master %>%
    dplyr::count(Population, IMPACT) %>%
    pivot_wider(names_from = IMPACT, values_from = n, values_fill = 0)

print(impact_by_pop)
write.csv(impact_by_pop, file.path(OUT, "impact_by_population.csv"), row.names = FALSE)

# =============================================================================
# 2D. Protein-coding variants only
# =============================================================================
cat("\n--- 2D. Protein-coding variants only ---\n")

coding <- master %>%
    filter(BIOTYPE == "protein_coding",
           IMPACT %in% c("HIGH", "MODERATE", "LOW")) %>%
    dplyr::select(Gene, Population, Drug, CHROM, POS, REF, ALT,
           Consequence, IMPACT, HGVSc, HGVSp, Existing_variation,
           AF_cohort, AN, ClinVar_CLNSIG, low_AN_flag, low_MQ_flag)

cat("Protein-coding non-MODIFIER variants:", nrow(coding), "\n")
print(coding %>% dplyr::select(Gene, Population, Consequence, IMPACT, HGVSp, AF_cohort, ClinVar_CLNSIG))
write.csv(coding, file.path(OUT, "coding_variants.csv"), row.names = FALSE)

# =============================================================================
# 2E. The single HIGH impact variant — inspect it
# =============================================================================
cat("\n--- 2E. HIGH impact variant(s) ---\n")

high_impact <- master[master$IMPACT == "HIGH", ]
cat("HIGH impact variants:", nrow(high_impact), "\n")
print(high_impact %>% dplyr::select(Gene, Population, CHROM, POS, REF, ALT,
                              Consequence, HGVSc, HGVSp, AF_cohort,
                              AN, ClinVar_CLNSIG, Existing_variation))
write.csv(high_impact, file.path(OUT, "high_impact_variants.csv"), row.names = FALSE)

# =============================================================================
# 2F. Plots
# =============================================================================
cat("\n--- 2F. Generating plots ---\n")

IMPACT_COLOURS <- c(
    HIGH     = "#d73027",
    MODERATE = "#fc8d59",
    LOW      = "#fee090",
    MODIFIER = "#e0f3f8"
)

# Plot 1: IMPACT stacked bar per gene
impact_long <- master %>%
    dplyr::count(Gene, Drug, IMPACT) %>%
    mutate(IMPACT = factor(IMPACT, levels = c("HIGH","MODERATE","LOW","MODIFIER")))

p1 <- ggplot(impact_long, aes(x = Gene, y = n, fill = IMPACT)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Drug, scales = "free_x") +
    scale_fill_manual(values = IMPACT_COLOURS) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Variant IMPACT Distribution per Gene",
         x = "Gene", y = "Number of variants", fill = "IMPACT")

ggsave(file.path(OUT, "plot_impact_by_gene.pdf"), p1, width = 10, height = 6)
ggsave(file.path(OUT, "plot_impact_by_gene.png"), p1, width = 10, height = 6, dpi = 150)

# Plot 2: Consequence breakdown
cons_plot <- master %>%
    mutate(IMPACT = factor(IMPACT, levels = c("HIGH","MODERATE","LOW","MODIFIER"))) %>%
    dplyr::count(Consequence_simple, IMPACT) %>%
    arrange(IMPACT, desc(n))

p2 <- ggplot(cons_plot, aes(x = reorder(Consequence_simple, n), y = n, fill = IMPACT)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = IMPACT_COLOURS) +
    theme_bw() +
    labs(title = "Consequence Type Distribution",
         x = "Consequence", y = "Count", fill = "IMPACT")

ggsave(file.path(OUT, "plot_consequence_distribution.pdf"), p2, width = 10, height = 6)
ggsave(file.path(OUT, "plot_consequence_distribution.png"), p2, width = 10, height = 6, dpi = 150)

cat("\nAnalysis 2 complete. Results saved to:", OUT, "\n")
