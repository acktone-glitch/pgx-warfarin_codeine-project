#!/bin/bash
#SBATCH --job-name=dbnsfp_download
#SBATCH --output=dbnsfp_download_%j.out
#SBATCH --error=dbnsfp_download_%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=15

set -euo pipefail

# Create directories
mkdir -p ~/reset/.vep/cache/dbNSFP
mkdir -p ~/reset/.vep/plugins

echo "Starting dbNSFP download..."
cd ~/reset/.vep/cache/dbNSFP

# Main data file with auto-retry
echo "Downloading main dbNSFP file..."
while ! curl --http1.1 -C - -O https://dist.genos.us/academic/e55b09/dbNSFP5.3.1a_grch38.gz; do
    echo "Download interrupted, retrying in 30 seconds..."
    sleep 30
done

# Index file with auto-retry
echo "Downloading index file..."
while ! curl --http1.1 -C - -O https://dist.genos.us/academic/e55b09/dbNSFP5.3.1a_grch38.gz.tbi; do
    echo "Download interrupted, retrying in 30 seconds..."
    sleep 30
done

# MD5 file with auto-retry
echo "Downloading md5 file..."
while ! curl --http1.1 -C - -O https://dist.genos.us/academic/e55b09/dbNSFP5.3.1a_grch38.gz.md5; do
    echo "Download interrupted, retrying in 30 seconds..."
    sleep 30
done

# Verify integrity
echo "Verifying download integrity..."
md5sum -c dbNSFP5.3.1a_grch38.gz.md5

# Download the VEP plugin
echo "Downloading dbNSFP VEP plugin..."
wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/115/dbNSFP.pm \
    -O ~/reset/.vep/plugins/dbNSFP.pm

if [[ -f ~/reset/.vep/plugins/dbNSFP.pm ]]; then
    echo "Plugin downloaded successfully."
else
    echo "ERROR: Plugin download failed." >&2
    exit 1
fi

echo "--------------------------------------"
echo "dbNSFP download and setup complete."
echo "Data: ~/.vep/cache/dbNSFP/dbNSFP5.3.1a_grch38.gz"
echo "Plugin: ~/.vep/plugins/dbNSFP.pm"
echo ""
echo "To use in VEP, add this to your VEP command:"
echo "--plugin dbNSFP,~/.vep/cache/dbNSFP/dbNSFP5.3.1a_grch38.gz,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,CADD_phred,REVEL_score,MutationTaster_pred,gnomAD_exomes_AF,LRT_score"
echo "--------------------------------------"
