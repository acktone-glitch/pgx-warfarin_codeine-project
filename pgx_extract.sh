#!/bin/bash
#SBATCH --job-name=pgx_extract_v3
#SBATCH --output=logs/pgx_extract_%j.log
#SBATCH --error=logs/pgx_extract_%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

set -euo pipefail

# ============================================================
#  ACTIVATE CONDA ENVIRONMENT (bcftools lives here)
# ============================================================
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null ||   source /etc/ace-data/home/rnagaba/miniconda3/etc/profile.d/conda.sh 2>/dev/null ||   source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
conda activate annotation 2>/dev/null || true

# ============================================================
#  CONFIGURATION
# ============================================================
BASE_DIR="/etc/ace-data/home/rnagaba/reset"
VARIANTS_DIR="${BASE_DIR}/variants"
POPULATIONS=("AFR" "EUR" "SAS" "EAS")
OUT_DIR="${BASE_DIR}/results/gene_vcfs"

# hg38 gene regions (no commas in coordinates)
declare -A GENE_REGION=(
  [CYP2C9]="chr10:96698415-96748853"
  [CYP2C19]="chr10:94762681-94855547"
  [CYP2D6]="chr22:42126499-42130865"
  [CYP4F2]="chr19:15878663-15916641"
  [VKORC1]="chr16:31096174-31111938"
  [GGCX]="chr2:85545546-85612953"
  [UGT2B7]="chr4:69048012-69118744"
  [OPRM1]="chr6:154039662-154137408"
  [DRD2]="chr11:113409605-113475553"
  [CALU]="chr7:128189484-128236610"
)

# ============================================================
#  SETUP
# ============================================================
mkdir -p "${OUT_DIR}" logs

# bcftools from PATH (no module needed)
BCFTOOLS=$(which bcftools 2>/dev/null || echo "")
if [[ -z "${BCFTOOLS}" ]]; then
  echo "ERROR: bcftools not found in PATH. Load module or provide full path."
  exit 1
fi
echo "bcftools: $("${BCFTOOLS}" --version | head -1)"
echo ""

# ============================================================
#  HELPER: detect compression and recompress to BGZF if needed
# ============================================================
# ============================================================
#  HELPER: check BGZF via magic bytes (bytes 12-13 = 66,67 = "BC")
#  Instant — reads only 16 bytes, no streaming needed
# ============================================================
is_bgzf() {
  local vcf="$1"
  python3 -c "
import sys
with open('${vcf}','rb') as f:
    h = f.read(16)
if len(h)>=16 and h[0]==0x1f and h[1]==0x8b and (h[3]&4) and h[12]==66 and h[13]==67:
    sys.exit(0)
sys.exit(1)
" 2>/dev/null
}

# ============================================================
#  HELPER: recompress to BGZF if needed, then index
# ============================================================
safe_index() {
  local vcf="$1"

  # Confirm BGZF
  if is_bgzf "${vcf}"; then
    echo "    ✓ BGZF confirmed"
  else
    echo "    FATAL: ${vcf} is not BGZF."
    echo "    Run manually: bgzip -@ 8 your_file.vcf  (rename .vcf first)"
    exit 1
  fi

  # Index with .csi (works for all chromosomes including chr22)
  if [[ ! -f "${vcf}.tbi" && ! -f "${vcf}.csi" ]]; then
    echo "    Indexing $(basename ${vcf}) with --csi..."
    "${BCFTOOLS}" index --csi --threads "${SLURM_CPUS_PER_TASK:-8}" "${vcf}"
    echo "    ✓ Indexed"
  else
    echo "    ✓ Index already exists"
  fi
}

# ============================================================
#  HELPER: chromosome inventory
# ============================================================
chrom_inventory() {
  local vcf="$1"
  "${BCFTOOLS}" index --stats "${vcf}" 2>/dev/null | \
    awk '{printf "    %-12s %s variants\n", $1, $3}' || \
    echo "    (could not read index stats)"
}


# ============================================================
#  CSQ FIELD POSITIONS from YOUR header (1-based after split on |)
#  Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|
#  BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|
#  Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|
#  STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|SOURCE|LoFtool|...
#  Then dbNSFP fields follow from position 25 onward
#
#  Key positions (1-based):
#   1  = Allele
#   2  = Consequence
#   3  = IMPACT
#   4  = SYMBOL         ← gene name filter
#   5  = Gene (Ensembl)
#   11 = HGVSc
#   12 = HGVSp
#   18 = Existing_variation (rsID etc)
#   25 = LoFtool
#
#  dbNSFP fields start at ~26. Key ones detected in your header:
#  Searching for exact positions of critical fields:
#   CADD_phred, REVEL_score, clinvar_clnsig, rs_dbSNP,
#   1000Gp3_AFR_AF, gnomAD4.1_joint_AFR_AF
#  These will be extracted by NAME using awk field-name lookup below.
# ============================================================

# AWK script that:
#  1. On first call, builds a field-name→position map from the CSQ FORMAT header
#  2. For each variant line, splits CSQ per transcript
#  3. Keeps the transcript matching the target gene (SYMBOL field)
#  4. Outputs a flat TSV row with named columns
# This avoids needing split-vep entirely.

AWK_EXTRACT='
BEGIN {
    OFS = "\t"
    # Fields we want to extract (in output order)
    split("SYMBOL,Consequence,IMPACT,Gene,BIOTYPE,EXON,INTRON,HGVSc,HGVSp," \
          "Existing_variation,LoFtool,CADD_phred,REVEL_score," \
          "AlphaMissense_score,AlphaMissense_pred," \
          "BayesDel_addAF_score,BayesDel_addAF_pred," \
          "SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred," \
          "MutationTaster_pred,clinvar_clnsig,clinvar_trait,clinvar_review," \
          "rs_dbSNP," \
          "1000Gp3_AFR_AF,1000Gp3_EAS_AF,1000Gp3_EUR_AF,1000Gp3_SAS_AF," \
          "gnomAD4.1_joint_AFR_AF,gnomAD4.1_joint_EAS_AF," \
          "gnomAD4.1_joint_NFE_AF,gnomAD4.1_joint_SAS_AF," \
          "ALFA_African_AF,AllofUs_AFR_AF,AllofUs_EAS_AF," \
          "AllofUs_EUR_AF,AllofUs_SAS_AF", \
          want_fields, ",")
    n_want = length(want_fields)
    header_printed = 0
    map_built = 0
}

# Build field position map from ##INFO=<ID=CSQ,...Format: ...> line
/^##INFO=<ID=CSQ/ {
    match($0, /Format: ([^"]+)/, arr)
    if (arr[1] != "") {
        n = split(arr[1], csq_names, "|")
        for (i = 1; i <= n; i++) {
            field_pos[csq_names[i]] = i
        }
        map_built = 1
    }
    next
}

# Skip other header lines
/^#/ { next }

# Process variant lines
{
    if (!map_built) {
        print "ERROR: CSQ format map not built — no ##INFO=<ID=CSQ line found" > "/dev/stderr"
        exit 1
    }

    chrom=$1; pos=$2; id=$3; ref=$4; alt=$5; qual=$6; filter=$7; info=$8

    # Extract AC and AN from INFO for cohort AF calculation
    ac = "."; an = "."; af_cohort = "."
    if (match(info, /\<AC=([^;]+)/, m)) ac = m[1]
    if (match(info, /\<AN=([^;]+)/, m)) an = m[1]
    if (ac != "." && an != "." && an+0 > 0) {
        af_cohort = sprintf("%.6f", ac / an)
    }

    # Find CSQ field in INFO
    csq_full = ""
    n_info = split(info, info_fields, ";")
    for (i = 1; i <= n_info; i++) {
        if (info_fields[i] ~ /^CSQ=/) {
            csq_full = substr(info_fields[i], 5)
            break
        }
    }
    if (csq_full == "") next

    # Split into per-transcript annotations
    n_csq = split(csq_full, csq_entries, ",")

    # Find best transcript for target gene:
    # Priority: CANONICAL=YES, else first match
    best = ""
    for (i = 1; i <= n_csq; i++) {
        n_f = split(csq_entries[i], f, "|")
        sym_pos = (4 in field_pos) ? field_pos["SYMBOL"] : 4
        if (f[sym_pos] == target_gene) {
            if (best == "") best = csq_entries[i]
            # Check if canonical
            can_pos = field_pos["VEP_canonical"]
            if (can_pos != "" && f[can_pos] == "YES") {
                best = csq_entries[i]; break
            }
        }
    }
    if (best == "") next  # variant not annotated for this gene

    # Print header once
    if (!header_printed) {
        printf "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tAC\tAN\tAF_cohort\tPopulation"
        for (i = 1; i <= n_want; i++) {
            printf "\t%s", want_fields[i]
        }
        printf "\n"
        header_printed = 1
    }

    # Split best transcript
    split(best, bf, "|")

    # Print variant fields
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", \
        chrom, pos, id, ref, alt, qual, filter, ac, an, af_cohort, population

    # Print CSQ fields by name lookup
    for (i = 1; i <= n_want; i++) {
        fname = want_fields[i]
        fpos  = field_pos[fname]
        val   = (fpos != "" && fpos <= length(bf)) ? bf[fpos] : "."
        printf "\t%s", val
    }
    printf "\n"
}
'

# ============================================================
#  SUMMARY TABLE
# ============================================================
SUMMARY="${OUT_DIR}/extraction_summary.tsv"
printf "%-6s\t%-10s\t%-10s\t%-8s\t%-10s\n" \
  "Pop" "Gene" "Total" "HIGH" "MODERATE" > "${SUMMARY}"

echo "================================================================"
echo "  PGx Gene VCF Extraction  v3"
echo "  $(date)"
echo "  Output: ${OUT_DIR}"
echo "================================================================"
echo ""

# ============================================================
#  PRE-FLIGHT: recompress to BGZF + index all population VCFs
#  Done once upfront so the main loop never stalls mid-job
# ============================================================
echo "── Pre-flight: BGZF check & indexing ──────────────────────"
for POP in "${POPULATIONS[@]}"; do
  VCF="${VARIANTS_DIR}/${POP}/${POP}_merged.vep.vcf.gz"
  if [[ ! -f "${VCF}" ]]; then
    echo "  ⚠  Missing: ${VCF} — ${POP} will be skipped"
    continue
  fi
  echo "  [${POP}] $(basename ${VCF})"
  safe_index "${VCF}"
  echo "  [${POP}] ready ✓"
done
echo ""

# ============================================================
#  MAIN LOOP: populations × genes
# ============================================================
for POP in "${POPULATIONS[@]}"; do

  VCF="${VARIANTS_DIR}/${POP}/${POP}_merged.vep.vcf.gz"
  POP_OUT="${OUT_DIR}/${POP}"
  mkdir -p "${POP_OUT}"

  echo "────────────────────────────────────────────────────────"
  echo "  Population: ${POP}"
  echo "  VCF: ${VCF}"

  if [[ ! -f "${VCF}" ]]; then
    echo "  ⚠  File not found — skipping ${POP}"; echo ""; continue
  fi

  # Chromosome inventory (index built in pre-flight)
  echo "  Chromosomes:"
  chrom_inventory "${VCF}"
  echo ""

  for GENE in "${!GENE_REGION[@]}"; do

    REGION="${GENE_REGION[$GENE]}"
    OUT_VCF="${POP_OUT}/${GENE}_${POP}.vcf.gz"
    OUT_TSV="${POP_OUT}/${GENE}_${POP}.tsv"

    # ── STEP 1: Extract gene VCF (3-tier fallback) ──────────

    # Tier 1: region + pipe-anchored symbol in CSQ
    "${BCFTOOLS}" view \
      --regions "${REGION}" \
      --threads "${SLURM_CPUS_PER_TASK:-4}" \
      "${VCF}" 2>/dev/null \
    | "${BCFTOOLS}" view \
        --include "INFO/CSQ~\"|${GENE}|\"" \
        --output-type z \
        --output "${OUT_VCF}" \
        --threads "${SLURM_CPUS_PER_TASK:-4}" \
        2>/dev/null \
    || true

    N_VARS=0
    [[ -f "${OUT_VCF}" ]] && \
      N_VARS=$("${BCFTOOLS}" view -H "${OUT_VCF}" 2>/dev/null | wc -l || echo 0)

    # Tier 2: symbol-only, no region (catches wider annotations)
    if [[ "${N_VARS}" -eq 0 ]]; then
      "${BCFTOOLS}" view \
        --include "INFO/CSQ~\"|${GENE}|\"" \
        --output-type z \
        --output "${OUT_VCF}" \
        --threads "${SLURM_CPUS_PER_TASK:-4}" \
        "${VCF}" 2>/dev/null || true
      [[ -f "${OUT_VCF}" ]] && \
        N_VARS=$("${BCFTOOLS}" view -H "${OUT_VCF}" 2>/dev/null | wc -l || echo 0)
    fi

    # Tier 3: loose match (no pipe anchors) — handles edge cases
    if [[ "${N_VARS}" -eq 0 ]]; then
      "${BCFTOOLS}" view \
        --regions "${REGION}" \
        --include "INFO/CSQ~\"${GENE}\"" \
        --output-type z \
        --output "${OUT_VCF}" \
        --threads "${SLURM_CPUS_PER_TASK:-4}" \
        "${VCF}" 2>/dev/null || true
      [[ -f "${OUT_VCF}" ]] && \
        N_VARS=$("${BCFTOOLS}" view -H "${OUT_VCF}" 2>/dev/null | wc -l || echo 0)
    fi

    if [[ "${N_VARS}" -eq 0 ]]; then
      echo "    ✗ ${GENE}: 0 variants in ${POP}"
      printf "%-6s\t%-10s\t%-10s\t%-8s\t%-10s\n" \
        "${POP}" "${GENE}" "0" "0" "0" >> "${SUMMARY}"
      continue
    fi

    # Index the gene VCF
    safe_index "${OUT_VCF}"

    # Count by IMPACT
    N_HIGH=$("${BCFTOOLS}" view -H "${OUT_VCF}" | \
      awk -F'\t' '{print $8}' | grep -oP '(?<=\|)HIGH(?=\|)' | wc -l || echo 0)
    N_MOD=$("${BCFTOOLS}" view -H "${OUT_VCF}" | \
      awk -F'\t' '{print $8}' | grep -oP '(?<=\|)MODERATE(?=\|)' | wc -l || echo 0)

    echo "    ✓ ${GENE} (${POP}): ${N_VARS} variants  [HIGH:${N_HIGH}  MOD:${N_MOD}]"

    # ── STEP 2: Flat TSV using awk CSQ parser ───────────────
    # Pass target gene and population as awk variables
    "${BCFTOOLS}" view "${OUT_VCF}" | \
      awk -v target_gene="${GENE}" -v population="${POP}" \
          "${AWK_EXTRACT}" \
    > "${OUT_TSV}"

    TSV_ROWS=$(wc -l < "${OUT_TSV}" || echo 0)
    echo "      TSV: ${TSV_ROWS} rows written → $(basename ${OUT_TSV})"

    printf "%-6s\t%-10s\t%-10s\t%-8s\t%-10s\n" \
      "${POP}" "${GENE}" "${N_VARS}" "${N_HIGH}" "${N_MOD}" >> "${SUMMARY}"

  done
  echo ""
done

# ============================================================
#  CYP2D6 DIAGNOSTIC
# ============================================================
echo "================================================================"
echo "  CYP2D6 Check"
echo "================================================================"
CYP2D6_TOTAL=0
for POP in "${POPULATIONS[@]}"; do
  F="${OUT_DIR}/${POP}/CYP2D6_${POP}.vcf.gz"
  N=0
  [[ -f "${F}" ]] && N=$("${BCFTOOLS}" view -H "${F}" 2>/dev/null | wc -l || echo 0)
  echo "  ${POP}: ${N} CYP2D6 variants"
  CYP2D6_TOTAL=$((CYP2D6_TOTAL + N))
done

if [[ "${CYP2D6_TOTAL}" -eq 0 ]]; then
  echo ""
  echo "  ⚠  CYP2D6 not found. Run these diagnostics:"
  echo ""
  echo "  # Is chr22 present in the VCF?"
  echo "  bcftools index --stats ${VARIANTS_DIR}/AFR/AFR_merged.vep.vcf.gz | grep chr22"
  echo ""
  echo "  # Any variants in the CYP2D6 window?"
  echo "  bcftools view -r chr22:42126499-42130865 \\"
  echo "      ${VARIANTS_DIR}/AFR/AFR_merged.vep.vcf.gz | grep -vc '^#'"
  echo ""
  echo "  # Any line with CYP2D6 in the INFO field?"
  echo "  bcftools view ${VARIANTS_DIR}/AFR/AFR_merged.vep.vcf.gz \\"
  echo "      | grep -v '^#' | grep -c 'CYP2D6' || echo 0"
fi

# ============================================================
#  FINAL SUMMARY
# ============================================================
echo ""
echo "================================================================"
echo "  SUMMARY — $(date)"
echo "================================================================"
column -t "${SUMMARY}"
echo ""
echo "  Gene VCFs : ${OUT_DIR}/<POP>/<GENE>_<POP>.vcf.gz"
echo "  Flat TSVs : ${OUT_DIR}/<POP>/<GENE>_<POP>.tsv"
echo "  Summary   : ${SUMMARY}"
echo "================================================================"
