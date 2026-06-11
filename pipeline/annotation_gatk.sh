#!/bin/bash
#SBATCH --job-name=vep_pgx
#SBATCH --output=/etc/ace-data/home/rnagaba/reset/logs/vep_%j.out
#SBATCH --error=/etc/ace-data/home/rnagaba/reset/logs/vep_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -eo pipefail

BASE="/etc/ace-data/home/rnagaba/reset"
CACHE_DIR="${BASE}/.vep/cache"
PLUGIN_DIR="${BASE}/.vep/plugins"
CLINVAR_VCF="${CACHE_DIR}/clinvar/clinvar.vcf.gz"
DBNSFP="${CACHE_DIR}/dbNSFP/dbNSFP5.3.1a_grch38.gz"
FASTA="${BASE}/reference/hg38.fa"
VEP="/etc/ace-data/home/rnagaba/.conda/envs/annotation/bin/vep"

POPS=(AFR EAS EUR SAS)

for POP in "${POPS[@]}"; do

    IN_VCF="${BASE}/results/gatk/${POP}/${POP}_pgx_filtered.vcf.gz"
    OUT_VCF="${BASE}/results/vep/${POP}/${POP}_pgx_filtered.vep.vcf.gz"
    STATS="${BASE}/results/vep/${POP}/${POP}_vep_stats.html"
    LOG="${BASE}/results/vep/${POP}/${POP}_vep.log"

    mkdir -p "${BASE}/results/vep/${POP}"

    if [[ ! -s "${IN_VCF}" ]]; then
        echo "❌ Missing VCF for ${POP} — skipping"
        continue
    fi

    echo "=== Annotating ${POP} ==="

    if ! ${VEP} \
        --cache \
        --offline \
        --dir_cache "${CACHE_DIR}" \
        --dir_plugins "${PLUGIN_DIR}" \
        --fasta "${FASTA}" \
        --species homo_sapiens \
        --assembly GRCh38 \
        --vcf \
        --input_file "${IN_VCF}" \
        --output_file "${OUT_VCF}" \
        --compress_output bgzip \
        --force_overwrite \
        --everything \
        --fork "${SLURM_CPUS_PER_TASK}" \
        --plugin LoFtool,"${PLUGIN_DIR}/LoFtool_scores.txt" \
         --plugin dbNSFP,"${DBNSFP}",\
SIFT_score,SIFT_pred,SIFT4G_score,SIFT4G_pred,\
Polyphen2_HDIV_score,Polyphen2_HDIV_pred,\
Polyphen2_HVAR_score,Polyphen2_HVAR_pred,\
MutationTaster_score,MutationTaster_pred,\
MetaSVM_score,MetaSVM_pred,\
MetaLR_score,MetaLR_pred,\
MetaRNN_score,MetaRNN_pred,\
CADD_phred,CADD_raw,\
REVEL_score,\
BayesDel_addAF_score,BayesDel_addAF_pred,\
BayesDel_noAF_score,BayesDel_noAF_pred,\
AlphaMissense_score,AlphaMissense_pred,\
GERP++_NR,GERP++_RS,\
phyloP100way_vertebrate,\
phastCons100way_vertebrate,\
1000Gp3_AF,1000Gp3_AFR_AF,1000Gp3_EAS_AF,1000Gp3_EUR_AF,1000Gp3_SAS_AF,\
gnomAD2.1.1_exomes_controls_AF,\
gnomAD2.1.1_exomes_controls_AFR_AF,\
gnomAD2.1.1_exomes_controls_EAS_AF,\
gnomAD2.1.1_exomes_controls_NFE_AF,\
gnomAD2.1.1_exomes_controls_SAS_AF,\
gnomAD4.1_joint_AF,\
gnomAD4.1_joint_AFR_AF,\
gnomAD4.1_joint_EAS_AF,\
gnomAD4.1_joint_NFE_AF,\
gnomAD4.1_joint_SAS_AF,\
clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review \
        --custom "${CLINVAR_VCF}",ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
        --stats_file "${STATS}" \
        2>&1 | tee "${LOG}"; then
            echo "❌ Error annotating ${POP} — see ${LOG}"
            continue
    fi

    echo "✅ Done ${POP} → ${OUT_VCF}"
done

echo "=== All populations complete ==="
