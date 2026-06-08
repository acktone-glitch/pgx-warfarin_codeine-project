#!/bin/bash
#SBATCH --job-name=download_fastq      # Name of the job
#SBATCH --output=download_%j.log       # Standard output and error log
#SBATCH --ntasks=2                     # Run a single task
#SBATCH --cpus-per-task=1              # Number of CPU cores per task
#SBATCH --mem=2G                       # Job memory request

# Download the paired-end fastq files
# Use -c to continue any partially downloaded files
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100031/SRR100031_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100031/SRR100031_2.fastq.gz
