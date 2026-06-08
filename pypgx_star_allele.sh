#!/bin/bash
#SBATCH --job-name=pypgx_pgx
#SBATCH --output=/etc/ace-data/home/rnagaba/reset/logs/pypgx_%A_%a.out
#SBATCH --error=/etc/ace-data/home/rnagaba/reset/logs/pypgx_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-23

# =============================================================================
#  PyPGx star allele calling — all 9 pharmacogenes × 4 populations
#  24 array tasks (6 genes × 4 pops = 24 combinations)
#
#  SUBMIT:
#    sbatch pypgx_star_alleles.sh
#
#  OUTPUT:
#    results/pypgx/<POP>/<GENE>/results.zip  — diplotype calls per sample
# =============================================================================

set -eo pipefail

# ── GENE × POPULATION MATRIX ─────────────────────────────────────────────────
# All 9 pharmacogenes — CYP2D6 included but treat with caution (short-read limits)
# GGCX, CALU, DRD2 not in PyPGx catalogue — PyPGx only covers P450/PGx-defined genes
GENES=(CYP2C9 CYP2D6 CYP4F2 VKORC1 UGT2B7 OPRM1)
POPS=(AFR EAS EUR SAS)

N_GENES=${#GENES[@]}   # 9
N_POPS=${#POPS[@]}     # 4
# Total tasks = 6 × 4 = 24  (array 0-23)

GENE_IDX=$(( SLURM_ARRAY_TASK_ID % N_GENES ))
POP_IDX=$(( SLURM_ARRAY_TASK_ID / N_GENES ))

GENE="${GENES[$GENE_IDX]}"
POP="${POPS[$POP_IDX]}"

BASE="/etc/ace-data/home/rnagaba/reset"
VCF="${BASE}/results/gatk/${POP}/${POP}_pgx_filtered.vcf.gz"
OUT="${BASE}/results/pypgx/${POP}/${GENE}"

echo "[INFO] Task:  ${SLURM_ARRAY_TASK_ID}"
echo "[INFO] Gene:  ${GENE}"
echo "[INFO] Pop:   ${POP}"
echo "[INFO] VCF:   ${VCF}"
echo "[INFO] Out:   ${OUT}"


# ── CHECK PYPGX ──────────────────────────────────────────────────────────────
if ! command -v pypgx &>/dev/null; then
    echo "[ERROR] pypgx not found — run: pip install pypgx"
    exit 1
fi
echo "[INFO] PyPGx: $(pypgx --version 2>&1 | head -1)"

# ── CHECK INPUT ──────────────────────────────────────────────────────────────
if [[ ! -f "${VCF}" ]]; then
    echo "[ERROR] VCF not found: ${VCF}"
    exit 1
fi

# ── SKIP IF ALREADY DONE ─────────────────────────────────────────────────────
if [[ -f "${OUT}/results.zip" ]]; then
    echo "[INFO] Already complete — skipping ${GENE}/${POP}"
    exit 0
fi

if [[ -d "${OUT}" ]]; then rm -rf "${OUT}"; fi

# ── RUN PYPGX ────────────────────────────────────────────────────────────────
pypgx run-ngs-pipeline \
    --variants "${VCF}" \
    --assembly GRCh38 \
    "${GENE}" \
    "${OUT}"

echo "[INFO] Done: ${GENE} in ${POP} → ${OUT}/results.zip"
