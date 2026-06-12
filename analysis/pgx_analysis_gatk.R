# =============================================================================
#  pgx_analysis_gatk.R
#  Full analysis suite for GATK-derived PGx master table
#  Analyses 01–08 matching the original bcftools analysis structure
#  but using MANE Select transcripts from the new GATK joint-called VCFs
#
#  USAGE: Open in RStudio, set ROOT below, then Source
#  INPUT: master_gatk_pass.csv (from parse_vep_gatk.R)
# =============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(stringr)

# ── PATHS ─────────────────────────────────────────────────────────────────────
ROOT    <- "D:/RITAH/SCHOOL/Beast Mode/analysis/analysis2"
IN_FILE <- file.path(ROOT, "results/master_gatk_pass.csv")

# Output directories — one per analysis
dirs <- file.path(ROOT, "results", paste0("0", 1:8, "_",
  c("summary","consequence","star_alleles","functional_impact",
    "clinvar","population_specific","loftool","cohort_correlation")))
for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
#  LOAD DATA
# =============================================================================
message("Loading: ", IN_FILE)
master <- read_csv(IN_FILE, show_col_types = FALSE)
message("Rows: ", nrow(master), "  Cols: ", ncol(master))

# Gene column name
gene_col <- if ("SYMBOL" %in% names(master)) "SYMBOL" else "Gene"

# Ensure Drug column exists
if (!"Drug" %in% names(master)) {
  master <- master %>% mutate(Drug = case_when(
    .data[[gene_col]] %in% c("CYP2C9","CYP2C19","VKORC1","CYP4F2","GGCX","CALU") ~ "warfarin",
    .data[[gene_col]] %in% c("CYP2D6","UGT2B7","OPRM1","DRD2")                   ~ "codeine",
    TRUE ~ "other"))
}

master <- master %>%
  rename(Gene = all_of(gene_col)) %>%
  mutate(
    AF_cohort = suppressWarnings(as.numeric(AF_cohort)),
    AN        = suppressWarnings(as.numeric(AN)),
    CADD_phred = suppressWarnings(as.numeric(CADD_phred)),
    low_AN_flag = !is.na(AN) & AN < 4,
    Consequence_simple = sapply(Consequence, function(x) {
      if (is.na(x)) return(NA)
      strsplit(x, "&")[[1]][1]
    })
  )

message("Genes: ", paste(sort(unique(master$Gene)), collapse=", "))
message("Pops:  ", paste(sort(unique(master$Population)), collapse=", "))

# =============================================================================
#  COLOUR PALETTES
# =============================================================================
POP_COLORS <- c(AFR="#D62728", EAS="#2CA02C", EUR="#1F77B4", SAS="#9467BD")
GENE_COLORS <- c(
  CALU="#E41A1C", CYP2C9="#377EB8", CYP2C19="#984EA3",
  CYP2D6="#4DAF4A", CYP4F2="#FF7F00", DRD2="#A65628",
  GGCX="#F781BF", OPRM1="#999999", UGT2B7="#FFFF33", VKORC1="#00CED1"
)
IMPACT_COLORS <- c(HIGH="#B2182B", MODERATE="#EF8A62", LOW="#FDD49E", MODIFIER="#D1E5F0")
DRUG_COLORS   <- c(warfarin="#2171B5", codeine="#CB181D", other="#AAAAAA")

theme_pgx <- function(base_size=12) {
  theme_bw(base_size=base_size) +
    theme(
      strip.background = element_rect(fill="grey92", colour="grey70"),
      strip.text       = element_text(face="bold"),
      panel.grid.major = element_line(colour="grey90"),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face="bold", size=base_size+1),
      plot.caption     = element_text(colour="grey50", size=base_size-3),
      axis.text        = element_text(colour="grey20")
    )
}

save_plot <- function(p, dir, name, w=10, h=6) {
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width=w, height=h)
  ggsave(file.path(dir, paste0(name, ".png")), p, width=w, height=h, dpi=150)
  message("  Saved: ", name)
}

# =============================================================================
#  ANALYSIS 01 — Summary statistics
# =============================================================================
message("\n=== Analysis 01: Summary ===")
OUT1 <- dirs[1]

summ <- master %>%
  group_by(Gene, Population, Drug) %>%
  summarise(
    n_variants   = n(),
    n_PASS       = sum(FILTER == "PASS", na.rm=TRUE),
    mean_AF      = round(mean(AF_cohort, na.rm=TRUE), 4),
    median_AF    = round(median(AF_cohort, na.rm=TRUE), 4),
    n_low_AN     = sum(low_AN_flag, na.rm=TRUE),
    .groups="drop"
  )
write_csv(summ, file.path(OUT1, "summary_by_gene_pop.csv"))

p1 <- ggplot(summ, aes(x=reorder(Gene, n_variants, sum),
                        y=n_variants, fill=Population)) +
  geom_bar(stat="identity", position="dodge", colour="white", linewidth=0.3) +
  coord_flip() +
  facet_wrap(~Drug, scales="free_y") +
  scale_fill_manual(values=POP_COLORS) +
  theme_pgx() +
  labs(title="Variant Count per Gene per Population (GATK PASS)",
       x="Gene", y="Number of variants",
       caption="MANE Select transcripts only")
save_plot(p1, OUT1, "plot_variant_counts", w=12, h=6)

message("  Summary table: ", nrow(summ), " rows")

# =============================================================================
#  ANALYSIS 02 — Consequence and IMPACT distribution
# =============================================================================
message("\n=== Analysis 02: Consequence / IMPACT ===")
OUT2 <- dirs[2]

# IMPACT stacked bar per gene
impact_long <- master %>%
  filter(!is.na(IMPACT)) %>%
  count(Gene, Drug, IMPACT) %>%
  mutate(IMPACT = factor(IMPACT, levels=c("HIGH","MODERATE","LOW","MODIFIER")))

p2a <- ggplot(impact_long, aes(x=Gene, y=n, fill=IMPACT)) +
  geom_bar(stat="identity", colour="white", linewidth=0.4) +
  coord_flip() +
  facet_wrap(~Drug, scales="free_x") +
  scale_fill_manual(values=IMPACT_COLORS) +
  theme_pgx() +
  labs(title="Variant IMPACT Distribution per Gene (GATK)",
       x="Gene", y="Number of variants")
save_plot(p2a, OUT2, "plot_impact_by_gene", w=11, h=6)

# Consequence type
cons_count <- master %>%
  filter(!is.na(Consequence_simple)) %>%
  count(Consequence_simple, IMPACT) %>%
  mutate(IMPACT = factor(IMPACT, levels=c("HIGH","MODERATE","LOW","MODIFIER"))) %>%
  arrange(desc(n))

p2b <- ggplot(cons_count,
              aes(x=reorder(Consequence_simple, n), y=n, fill=IMPACT)) +
  geom_bar(stat="identity", colour="white", linewidth=0.4) +
  coord_flip() +
  scale_fill_manual(values=IMPACT_COLORS) +
  theme_pgx() +
  labs(title="Consequence Type Distribution (GATK)",
       x="Consequence", y="Count")
save_plot(p2b, OUT2, "plot_consequence_distribution", w=10, h=6)

# IMPACT per population
p2c <- master %>%
  filter(!is.na(IMPACT)) %>%
  count(Population, IMPACT) %>%
  mutate(IMPACT = factor(IMPACT, levels=c("HIGH","MODERATE","LOW","MODIFIER"))) %>%
  ggplot(aes(x=Population, y=n, fill=IMPACT)) +
  geom_bar(stat="identity", colour="white", linewidth=0.4) +
  scale_fill_manual(values=IMPACT_COLORS) +
  theme_pgx() +
  labs(title="IMPACT Distribution per Population (GATK)",
       x="Population", y="Number of variants")
save_plot(p2c, OUT2, "plot_impact_by_population", w=8, h=5)

write_csv(cons_count, file.path(OUT2, "consequence_counts.csv"))

# =============================================================================
#  ANALYSIS 03 — Star allele identification (rsID-based)
# =============================================================================
message("\n=== Analysis 03: Star alleles ===")
OUT3 <- dirs[3]

star_catalogue <- data.frame(
  Gene           = c("CYP2C9","CYP2C9","CYP2C9","CYP2C9","CYP2C9","CYP2C9",
                     "CYP2C19","CYP2C19","CYP2C19","CYP4F2","VKORC1","UGT2B7","OPRM1"),
  Star_allele    = c("*2","*3","*5","*6","*8","*11",
                     "*2","*3","*17","*3","H1","*2","118A>G"),
  rsID           = c("rs1799853","rs1057910","rs28371686","rs9332131","rs7900194","rs28371685",
                     "rs4244285","rs4986893","rs12248560","rs2108622","rs9923231","rs7439366","rs1799971"),
  HGVSp          = c("p.Arg144Cys","p.Ile359Leu","p.Asp360Glu","p.Lys395del","p.Arg150His","p.Arg335Trp",
                     "splice defect","p.Trp212Ter","promoter -806C>T","p.Val433Met","c.-1639G>A","p.His268Tyr","p.Asn40Asp"),
  Expected_effect= c("Reduced function","No function","No function","No function","Reduced function","Reduced function",
                     "No function","No function","Increased function","Reduced function","Warfarin sensitive","Normal function","Reduced expression"),
  stringsAsFactors = FALSE
)

# Match by rsID (Existing_variation column) or HGVSp
ev_col <- if ("Existing_variation" %in% names(master)) "Existing_variation" else
          if ("ID" %in% names(master)) "ID" else NA

hgvsp_col <- if ("HGVSp" %in% names(master)) "HGVSp" else NA

star_found <- do.call(rbind, lapply(seq_len(nrow(star_catalogue)), function(i) {
  row <- star_catalogue[i,]
  hits <- master %>%
    filter(Gene == row$Gene) %>%
    filter(
      (!is.na(ev_col) & str_detect(coalesce(.data[[ev_col]], ""), fixed(row$rsID))) |
      (!is.na(hgvsp_col) & str_detect(coalesce(.data[[hgvsp_col]], ""), fixed(row$HGVSp)))
    )
  if (nrow(hits) == 0) return(NULL)
  hits %>%
    group_by(Population) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(Star_allele   = row$Star_allele,
           rsID_target   = row$rsID,
           Expected_effect = row$Expected_effect)
}))

write_csv(star_catalogue, file.path(OUT3, "star_allele_catalogue.csv"))

if (!is.null(star_found) && nrow(star_found) > 0) {
  write_csv(star_found, file.path(OUT3, "star_alleles_found.csv"))
  message("  Found ", nrow(star_found), " star allele hits across populations")

  p3 <- star_found %>%
    filter(!is.na(AF_cohort)) %>%
    mutate(Label = paste(Gene, Star_allele)) %>%
    ggplot(aes(x=Population, y=AF_cohort, fill=Population,
               label=round(AF_cohort, 3))) +
    geom_bar(stat="identity", colour="white", linewidth=0.4) +
    geom_text(vjust=-0.3, size=2.8, colour="grey20") +
    facet_wrap(~Label, scales="free_y") +
    scale_fill_manual(values=POP_COLORS) +
    theme_pgx() +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.position="none") +
    labs(title="PGx Star Allele Frequencies by Population (GATK)",
         x="Population", y="Allele Frequency (cohort)",
         caption="Matched by rsID or HGVSp notation")
  save_plot(p3, OUT3, "plot_star_allele_frequencies", w=14, h=10)
} else {
  message("  No star alleles matched — check Existing_variation column name")
  message("  Available columns: ", paste(names(master)[1:20], collapse=", "))
}

# =============================================================================
#  ANALYSIS 04 — Functional impact / CADD / pathogenicity
# =============================================================================
message("\n=== Analysis 04: Functional impact ===")
OUT4 <- dirs[4]

missense <- master %>%
  filter(IMPACT %in% c("MODERATE","HIGH")) %>%
  mutate(
    CADD_phred = suppressWarnings(as.numeric(CADD_phred)),
    SIFT_num   = suppressWarnings(as.numeric(str_extract(
                   coalesce(SIFT_pred, SIFT_score, ""), "^[0-9.]+"
                 ))),
    PolyPhen_num = suppressWarnings(as.numeric(str_extract(
                   coalesce(Polyphen2_HDIV_score, PolyPhen, ""), "^[0-9.]+"
                 ))),
    REVEL_num  = suppressWarnings(as.numeric(REVEL_score)),
    CADD_flag  = !is.na(CADD_phred) & CADD_phred >= 20,
    SIFT_flag  = !is.na(SIFT_num)   & SIFT_num   <= 0.05,
    PP2_flag   = !is.na(PolyPhen_num) & PolyPhen_num >= 0.85,
    REVEL_flag = !is.na(REVEL_num)  & REVEL_num  >= 0.5,
    n_tools_flagged = rowSums(cbind(CADD_flag, SIFT_flag, PP2_flag, REVEL_flag),
                              na.rm=TRUE),
    Pathogenicity_class = case_when(
      IMPACT == "HIGH"               ~ "HIGH impact (no scores)",
      n_tools_flagged >= 3           ~ "Likely damaging",
      n_tools_flagged == 2           ~ "Possibly damaging",
      n_tools_flagged == 1           ~ "Uncertain",
      TRUE                           ~ "Likely benign / no data"
    )
  )

write_csv(missense, file.path(OUT4, "missense_functional_scores.csv"))
message("  MODERATE/HIGH variants: ", nrow(missense))

PATHO_COLORS <- c(
  "Likely damaging"        = "#B2182B",
  "Possibly damaging"      = "#EF8A62",
  "Uncertain"              = "#FDD49E",
  "Likely benign / no data"= "#92C5DE",
  "HIGH impact (no scores)"= "#4D0000"
)

# CADD distribution
if (any(!is.na(missense$CADD_phred))) {
  p4a <- missense %>%
    filter(!is.na(CADD_phred)) %>%
    ggplot(aes(x=CADD_phred, fill=Gene)) +
    geom_histogram(bins=20, colour="white", linewidth=0.3) +
    geom_vline(xintercept=20, linetype="dashed",
               colour="#B2182B", linewidth=0.8) +
    annotate("text", x=21, y=Inf, label="CADD≥20\n(top 1%)",
             hjust=0, vjust=1.5, colour="#B2182B", size=3.2) +
    facet_wrap(~Drug) +
    scale_fill_manual(values=GENE_COLORS, na.value="grey70") +
    theme_pgx() +
    labs(title="CADD Score Distribution — Moderate/High Impact Variants (GATK)",
         x="CADD phred score", y="Count",
         caption="Red dashed line = CADD 20 threshold (top 1% most deleterious)")
  save_plot(p4a, OUT4, "plot_cadd_distribution", w=10, h=5)
}

# Pathogenicity classification
p4b <- ggplot(missense,
              aes(x=reorder(Gene, Gene, function(x) length(x)),
                  fill=Pathogenicity_class)) +
  geom_bar(colour="white", linewidth=0.4) +
  coord_flip() +
  scale_fill_manual(values=PATHO_COLORS) +
  theme_pgx() +
  labs(title="Missense Variant Pathogenicity Classification per Gene (GATK)",
       x="Gene", y="Number of variants")
save_plot(p4b, OUT4, "plot_pathogenicity_class", w=10, h=5)

# =============================================================================
#  ANALYSIS 05 — ClinVar drug response
# =============================================================================
message("\n=== Analysis 05: ClinVar ===")
OUT5 <- dirs[5]

# Identify ClinVar column
cv_col <- names(master)[str_detect(names(master), regex("clin.*sig|CLIN_SIG|ClinVar", ignore_case=TRUE))][1]

if (!is.na(cv_col)) {
  message("  ClinVar column: ", cv_col)
  clinvar_all <- master %>%
    filter(!is.na(.data[[cv_col]])) %>%
    rename(ClinVar_CLNSIG = all_of(cv_col))

  drug_response <- clinvar_all %>%
    filter(str_detect(ClinVar_CLNSIG, "drug_response"))

  write_csv(clinvar_all,    file.path(OUT5, "clinvar_all_variants.csv"))
  write_csv(drug_response,  file.path(OUT5, "drug_response_variants.csv"))
  message("  ClinVar entries: ", nrow(clinvar_all),
          "  Drug response: ", nrow(drug_response))

  CLINVAR_COLORS <- c(
    "drug_response"         = "#D62728",
    "Benign/Likely_benign"  = "#2CA02C",
    "Uncertain_significance"= "#FF7F0E",
    "Pathogenic"            = "#9467BD",
    "Other"                 = "#AAAAAA"
  )

  p5a <- clinvar_all %>%
    mutate(CLNSIG_simple = case_when(
      str_detect(ClinVar_CLNSIG, "drug_response") ~ "drug_response",
      str_detect(ClinVar_CLNSIG, "Pathogenic")    ~ "Pathogenic",
      str_detect(ClinVar_CLNSIG, "Benign")        ~ "Benign/Likely_benign",
      str_detect(ClinVar_CLNSIG, "Uncertain")     ~ "Uncertain_significance",
      TRUE ~ "Other"
    )) %>%
    count(Gene, CLNSIG_simple) %>%
    ggplot(aes(x=reorder(Gene, n, sum), y=n, fill=CLNSIG_simple)) +
    geom_bar(stat="identity", colour="white", linewidth=0.4) +
    coord_flip() +
    scale_fill_manual(values=CLINVAR_COLORS, name="ClinVar significance") +
    theme_pgx() +
    labs(title="ClinVar Clinical Significance per Gene (GATK)",
         x="Gene", y="Number of variants")
  save_plot(p5a, OUT5, "plot_clinvar_by_gene", w=9, h=6)

  if (nrow(drug_response) > 0) {
    p5b <- drug_response %>%
      filter(!is.na(AF_cohort)) %>%
      ggplot(aes(x=Population, y=AF_cohort, fill=Population)) +
      geom_boxplot(alpha=0.75, outlier.shape=NA) +
      geom_jitter(aes(colour=Population), width=0.18, size=1.6, alpha=0.7) +
      facet_wrap(~Gene, scales="free_y") +
      scale_fill_manual(values=POP_COLORS) +
      scale_colour_manual(values=POP_COLORS) +
      theme_pgx() + theme(legend.position="none") +
      labs(title="Drug Response Variant AF by Population (GATK)",
           x="Population", y="Allele Frequency (cohort)")
    save_plot(p5b, OUT5, "plot_drug_response_af", w=12, h=8)
  }
} else {
  message("  WARNING: No ClinVar column found. Available: ",
          paste(names(master)[1:30], collapse=", "))
}

# =============================================================================
#  ANALYSIS 06 — Population-specific variants
# =============================================================================
message("\n=== Analysis 06: Population-specific ===")
OUT6 <- dirs[6]

# Identify gnomAD / 1000G population AF columns
afr_col <- names(master)[str_detect(names(master), "AFR_AF") &
                          str_detect(names(master), "1000G|1000Gp3|gnomAD4")][1]
eas_col <- names(master)[str_detect(names(master), "EAS_AF") &
                          str_detect(names(master), "1000G|1000Gp3|gnomAD4")][1]
eur_col <- names(master)[str_detect(names(master), "(?:NFE|EUR)_AF") &
                          str_detect(names(master), "1000G|1000Gp3|gnomAD4")][1]
sas_col <- names(master)[str_detect(names(master), "SAS_AF") &
                          str_detect(names(master), "1000G|1000Gp3|gnomAD4")][1]

message("  Reference AF columns: AFR=", afr_col, " EAS=", eas_col,
        " EUR=", eur_col, " SAS=", sas_col)

# AFR-enriched variants
if (!is.na(afr_col) && !is.na(eas_col) && !is.na(eur_col) && !is.na(sas_col)) {
  pop_specific <- master %>%
    mutate(
      ref_AFR = suppressWarnings(as.numeric(.data[[afr_col]])),
      ref_EAS = suppressWarnings(as.numeric(.data[[eas_col]])),
      ref_EUR = suppressWarnings(as.numeric(.data[[eur_col]])),
      ref_SAS = suppressWarnings(as.numeric(.data[[sas_col]]))
    ) %>%
    mutate(
      AFR_enriched = !is.na(ref_AFR) & ref_AFR > 0.05 &
                     (is.na(ref_EAS) | ref_AFR > 2*ref_EAS) &
                     (is.na(ref_EUR) | ref_AFR > 2*ref_EUR) &
                     (is.na(ref_SAS) | ref_AFR > 2*ref_SAS),
      ref_rare = (!is.na(AF_cohort) & AF_cohort > 0.01) &
                 (is.na(ref_AFR) | ref_AFR < 0.001)
    )

  write_csv(pop_specific %>% filter(AFR_enriched),
            file.path(OUT6, "afr_enriched_variants.csv"))
  write_csv(pop_specific %>% filter(ref_rare),
            file.path(OUT6, "cohort_common_ref_rare.csv"))

  message("  AFR-enriched: ", sum(pop_specific$AFR_enriched, na.rm=TRUE))
  message("  Cohort-common/ref-rare: ", sum(pop_specific$ref_rare, na.rm=TRUE))

  # AFR vs EUR scatter
  scatter_dat <- pop_specific %>%
    filter(!is.na(ref_AFR), !is.na(ref_EUR))

  if (nrow(scatter_dat) > 0) {
    p6a <- ggplot(scatter_dat, aes(x=ref_EUR, y=ref_AFR, colour=Gene)) +
      geom_point(alpha=0.65, size=2) +
      geom_abline(slope=1, intercept=0,
                  linetype="dashed", colour="grey40", linewidth=0.7) +
      scale_colour_manual(values=GENE_COLORS, na.value="grey70") +
      theme_pgx() +
      labs(title="AFR vs EUR Reference Allele Frequency (GATK variants)",
           x="EUR allele frequency", y="AFR allele frequency",
           caption="Points above diagonal = higher in AFR")
    save_plot(p6a, OUT6, "plot_afr_vs_eur_af", w=8, h=6)
  }
}

# Cohort AF vs gnomAD4 scatter
gnom_col <- names(master)[str_detect(names(master), "gnomAD4.1_joint_AF$")][1]
if (!is.na(gnom_col)) {
  cohort_vs_ref <- master %>%
    filter(!is.na(AF_cohort)) %>%
    mutate(ref_gnomad = suppressWarnings(as.numeric(.data[[gnom_col]]))) %>%
    filter(!is.na(ref_gnomad))

  p6b <- ggplot(cohort_vs_ref,
                aes(x=ref_gnomad, y=AF_cohort, colour=Population)) +
    geom_point(alpha=0.5, size=1.6) +
    geom_abline(slope=1, intercept=0,
                linetype="dashed", colour="grey40", linewidth=0.7) +
    facet_wrap(~Population) +
    scale_colour_manual(values=POP_COLORS) +
    theme_pgx() + theme(legend.position="none") +
    labs(title="Cohort AF vs gnomAD 4.1 AF (GATK variants)",
         x="gnomAD 4.1 joint AF", y="Cohort AF",
         caption="Dashed = perfect agreement")
  save_plot(p6b, OUT6, "plot_cohort_vs_gnomad4", w=10, h=8)
}

# =============================================================================
#  ANALYSIS 07 — LoFtool gene sensitivity
# =============================================================================
message("\n=== Analysis 07: LoFtool ===")
OUT7 <- dirs[7]

loftool_col <- names(master)[str_detect(names(master), regex("loftool", ignore_case=TRUE))][1]

if (!is.na(loftool_col)) {
  loftool_summary <- master %>%
    mutate(LoFtool_score = suppressWarnings(as.numeric(.data[[loftool_col]]))) %>%
    filter(!is.na(LoFtool_score)) %>%
    group_by(Gene, Drug) %>%
    summarise(LoFtool_score = first(LoFtool_score), .groups="drop")

  write_csv(loftool_summary, file.path(OUT7, "loftool_scores.csv"))

  p7 <- ggplot(loftool_summary,
               aes(x=reorder(Gene, LoFtool_score), y=LoFtool_score, fill=Drug)) +
    geom_bar(stat="identity", colour="white", linewidth=0.4) +
    geom_hline(yintercept=0.5, linetype="dashed",
               colour="grey40", linewidth=0.7) +
    coord_flip() +
    scale_fill_manual(values=DRUG_COLORS) +
    theme_pgx() +
    labs(title="LoFtool Score per Gene",
         subtitle="Lower = more intolerant to loss-of-function variants",
         x="Gene", y="LoFtool score",
         caption="Dashed line = 0.5 midpoint")
  save_plot(p7, OUT7, "plot_loftool_by_gene", w=8, h=5)
} else {
  message("  WARNING: No LoFtool column found")
}

# =============================================================================
#  ANALYSIS 08 — Cohort AF concordance with reference populations
# =============================================================================
message("\n=== Analysis 08: Cohort concordance ===")
OUT8 <- dirs[8]

ref_map <- list(
  AFR = names(master)[str_detect(names(master), "AFR_AF") &
                       str_detect(names(master), "1000Gp3|gnomAD4")][1],
  EAS = names(master)[str_detect(names(master), "EAS_AF") &
                       str_detect(names(master), "1000Gp3|gnomAD4")][1],
  EUR = names(master)[str_detect(names(master), "(?:NFE|EUR)_AF") &
                       str_detect(names(master), "1000Gp3|gnomAD4")][1],
  SAS = names(master)[str_detect(names(master), "SAS_AF") &
                       str_detect(names(master), "1000Gp3|gnomAD4")][1]
)

corr_results <- do.call(rbind, lapply(names(ref_map), function(pop) {
  ref_col <- ref_map[[pop]]
  if (is.na(ref_col)) return(NULL)

  dat <- master %>%
    filter(Population == pop, !is.na(AF_cohort)) %>%
    mutate(ref_AF = suppressWarnings(as.numeric(.data[[ref_col]]))) %>%
    filter(!is.na(ref_AF))

  if (nrow(dat) < 5) return(NULL)

  ct <- cor.test(dat$AF_cohort, dat$ref_AF, method="pearson")
  data.frame(Population=pop, r=round(ct$estimate, 3),
         p_value=ct$p.value, n=nrow(dat), ref_column=ref_col)
}))

if (!is.null(corr_results) && nrow(corr_results) > 0) {
  write_csv(corr_results, file.path(OUT8, "cohort_reference_correlations.csv"))
  message("  Correlation results:")
  print(corr_results)

  # Scatter per population with correlation annotation
  scatter_all <- do.call(rbind, lapply(names(ref_map), function(pop) {
    ref_col <- ref_map[[pop]]
    if (is.na(ref_col)) return(NULL)
    master %>%
      filter(Population == pop, !is.na(AF_cohort)) %>%
      mutate(ref_AF = suppressWarnings(as.numeric(.data[[ref_col]])),
             Population = pop) %>%
      filter(!is.na(ref_AF)) %>%
      select(Population, AF_cohort, ref_AF, Gene)
  }))

  r_labels <- corr_results %>%
    mutate(label = paste0("r = ", r, "\n(n=", n, ")"))

  p8 <- ggplot(scatter_all, aes(x=ref_AF, y=AF_cohort, colour=Population)) +
    geom_point(alpha=0.5, size=1.6) +
    geom_abline(slope=1, intercept=0,
                linetype="dashed", colour="grey40", linewidth=0.7) +
    geom_text(data=r_labels,
              aes(x=0.7, y=0.1, label=label),
              inherit.aes=FALSE, size=3, colour="grey30") +
    facet_wrap(~Population) +
    scale_colour_manual(values=POP_COLORS) +
    theme_pgx() + theme(legend.position="none") +
    labs(title="Cohort AF vs Reference Population AF (GATK variants)",
         x="Reference AF (1000G or gnomAD4)",
         y="Cohort AF",
         caption="Dashed = perfect agreement")
  save_plot(p8, OUT8, "plot_cohort_vs_reference", w=10, h=8)
}

# =============================================================================
#  FINAL SUMMARY
# =============================================================================
message("\n", paste(rep("=", 60), collapse=""))
message("ALL ANALYSES COMPLETE")
message(paste(rep("=", 60), collapse=""))
message("Output directories:")
for (d in dirs) message("  ", basename(d), "/")
message("")
message("Key counts:")
message("  Total PASS variants (MANE Select): ", nrow(master))
message("  Unique genes:      ", length(unique(master$Gene)))
message("  Populations:       ", paste(sort(unique(master$Population)), collapse=", "))
if (exists("star_found") && !is.null(star_found))
  message("  Star alleles found: ", nrow(star_found))
if (exists("clinvar_all") && !is.null(clinvar_all))
  message("  ClinVar entries:   ", nrow(clinvar_all))
if (exists("drug_response") && !is.null(drug_response))
  message("  Drug response:     ", nrow(drug_response))
