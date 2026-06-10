# Pharmacogenomics of Warfarin and Codeine Metabolism Across Global Populations

**Internship project | Infectious Disease Institute, Makerere University**  
**SHEDS Programme / African Centre of Excellence in Bioinformatics (ACE)**  
**Supervisor: Kakembo Elishmah Fredrick**

---

## Overview

This repository contains the full bioinformatics pipeline, R analysis scripts, and summary results for a pharmacogenomics study characterising genetic variation in warfarin and codeine metabolism genes across four ancestral populations (AFR, EAS, EUR, SAS) using whole-exome sequencing data from the 1000 Genomes Project.

**Drug pathways studied:**
- **Warfarin:** CYP2C9, CYP2C19, CYP4F2, VKORC1, GGCX, CALU
- **Codeine:** CYP2D6, UGT2B7, OPRM1, DRD2

**Key findings:**
- 3,657 PASS variants on canonical (MANE Select) transcripts across 4 populations
- CYP2C19 reduced-function metabolisers prevalent in SAS (45.5%) and AFR (41.7%)
- CYP4F2*3 at clinically actionable frequencies in SAS (41%) and EUR (30%)
- Cohort allele frequencies highly concordant with 1000G reference (r ≥ 0.942)
- 5 of 13 queried star alleles recovered from VEP annotation

---

## Repository Structure

```
├── pipeline/           SLURM bash scripts for HPC variant calling & annotation
├── analysis/           R scripts for statistical analysis and visualisation
├── config/             PGx target intervals (BED format)
├── results/            Summary tables (CSV) — no raw VCFs
└── docs/               Results report and figures
```

---

## Pipeline

### Requirements
| Tool | Version | Purpose |
| --- | --- | --- |
| GATK | 4.6.2.0 | Joint variant calling |
| VEP | 115 | Variant annotation |
| PyPGx | 0.26.0 | Star allele calling |
| bcftools | 1.x | VCF manipulation |
| samtools | 1.x | BAM processing |
| fastp | 0.23+ | Read QC and trimming |
| R | 4.5+ | Statistical analysis |
### Steps

```
00_fastp_qc.sh               FASTQ quality control + adapter trimming
01_bwa_mem2_align.sh         Alignment to GRCh38 using BWA‑MEM2
02_markduplicates.sh         Mark PCR duplicates (GATK MarkDuplicates)
03_bqsr.sh                   Base Quality Score Recalibration (BQSR)

04_gatk_haplotypecaller.sh   Per-sample GVCF generation (SLURM array)
05_gatk_genomicsdb.sh        GenomicsDBImport across PGx intervals
06_gatk_genotype.sh          Joint genotyping + hard filtering

07_vep_annotation.sh         VEP annotation with dbNSFP, ClinVar, LoFtool
08_pypgx_star_alleles.sh     PyPGx haplotype-resolved star allele calling
09_compare_gatk_bcftools.sh  GATK vs bcftools variant caller concordance

```

### HPC Configuration

- Cluster: `kla-ac-hpc-02` (SLURM scheduler)
- Reference: GRCh38 (hg38)
- Conda environment: `annotation`
- Base path: `/etc/ace-data/home/rnagaba/reset/`

---

## Analysis

All R scripts are in `analysis/` and run in RStudio on Windows.

```
parse_vep_gatk.R       Parse VEP-annotated VCF files into master table
pgx_analysis_gatk.R    8-part analysis: consequences, star alleles, ClinVar,
                        CADD, population-specific variants, LoFtool, concordance
pypgx_analysis.R       PyPGx diplotype frequencies, phenotype distributions,
                        activity scores, drug pathway profiles
```

---

## Key Results

| Gene | Key variant | AFR AF | EAS AF | EUR AF | SAS AF |
|------|------------|--------|--------|--------|--------|
| CYP2C19 | *2 (rs4244285) | 0.221 | 0.056 | 0.146 | 0.318 |
| CYP2C19 | *3 (rs4986893) | — | 0.050 | — | — |
| CYP4F2 | *3 (rs2108622) | 0.076 | 0.182 | 0.300 | 0.409 |
| OPRM1 | 118A>G (rs1799971) | 0.013 | 0.333 | 0.167 | 0.400 |
| UGT2B7 | *2 (rs7439366) | 0.821 | 0.611 | 0.479 | 0.545 |

**CYP2C19 reduced-function phenotype burden (PyPGx):**

| Population | Poor Metabolizer | Intermediate Metabolizer | Combined |
|-----------|-----------------|-------------------------|---------|
| AFR | 8.3% | 33.3% | **41.7%** |
| EAS | 16.7% | 0% | **16.7%** |
| EUR | 4.0% | 24.0% | **28.0%** |
| SAS | 9.1% | 36.4% | **45.5%** |

---

## Data Availability

Raw WES data: 1000 Genomes Project (https://www.internationalgenome.org/)  
Reference genome: GRCh38/hg38  
VEP cache: Ensembl v115  
dbNSFP: v5.3.1a

Raw VCF files are not included in this repository (size limits).  
Summary tables are provided in `results/`.

---

## Citation

If you use this pipeline or results, please cite:

> Nagaba Ritah Acktone (2025). Pharmacogenomics of warfarin and codeine metabolism  
> across global populations using whole-exome sequencing data.  
> Internship Report, Infectious Disease Institute, Makerere University.

---

## Contact

**Nagaba Ritah Acktone**  
ritahtoni@gmail.com
Human Genomics Intern, IDI/ACE/SHEDS Programme  
Makerere University, Kampala, Uganda
