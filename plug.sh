#!/bin/bash
#SBATCH --job-name=vep_install
#SBATCH --output=vep_install_%j.out
#SBATCH --error=vep_install_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
set -euo pipefail

# Load Perl if needed (depends on your cluster)
# module load perl

# Move into the VEP installation directory
cd ensembl-vep

# Create cache + plugin directories
mkdir -p ~/reset/.vep/cache
mkdir -p ~/reset/.vep/plugins

# Run the VEP installer
perl INSTALL.pl \
  --NO_UPDATE \
  --AUTO cfp \
  --SPECIES homo_sapiens \
  --ASSEMBLY GRCh38 \
  --CACHEDIR ~/reset/.vep/cache \
  --PLUGINSDIR ~/reset/.vep/plugins \
  --CACHE_VERSION 115 \
  --PLUGINS LoFtool

echo "VEP cache + plugin installation complete."
echo "Remember: ClinVar and gnomAD must be downloaded manually and indexed."
