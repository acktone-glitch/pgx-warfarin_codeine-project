#!/bin/bash
#SBATCH --job-name=index_hg38
#SBATCH --output=logs/index_hg38_%A.out
#SBATCH --error=logs/index_hg38_%A.err
#SBATCH --cpus-per-task=48
#SBATCH --mem=32G

source ~/.bashrc
conda activate alignment

cd reference

echo "Indexing hg38 with BWA..."
bwa index hg38

echo "Indexing hg38 with samtools..."
samtools faidx hg38

echo "Creating sequence dictionary..."
samtools dict hg38 -o hg38.dict

echo "Reference indexing complete."
