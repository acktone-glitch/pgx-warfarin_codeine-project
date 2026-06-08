#!/bin/bash
#SBATCH --job-name=flagstat_all
#SBATCH --output=logs/flagstat_%A.out
#SBATCH --error=logs/flagstat_%A.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G


# Starting directory (reset/)
BASE_DIR="/etc/ace-data/home/rnagaba/reset"
DATA_DIR="$BASE_DIR/data"

# Output file
OUT="$BASE_DIR/flagstat_summary.tsv"
mkdir -p $BASE_DIR/logs

echo -e "Sample\tTotal_Reads\tMapped_Reads\tPercent_Mapped" > $OUT

echo "Searching for BAM files under: $DATA_DIR"

# Loop through all BAM files under all superpops → pops → samples → aligned/
find $DATA_DIR -type f -name "*.bam" | while read bam; do
    sample=$(basename "$bam" .bam)

    echo "Processing: $sample"

    stats=$(samtools flagstat "$bam")

    total=$(echo "$stats" | grep "in total" | awk '{print $1}')
    mapped=$(echo "$stats" | grep "mapped (" | awk '{print $1}')
    percent=$(echo "$stats" | grep "mapped (" | awk -F'[()%]' '{print $2}')

    echo -e "${sample}\t${total}\t${mapped}\t${percent}" >> $OUT
done

echo "Flagstat summary written to $OUT"

# Compute averages
echo -e "\nAVERAGES:" >> $OUT
awk 'NR>1 && $2>0 {total+=$2; mapped+=$3; pct+=$4; n++}
     END {
         print "Avg_Total_Reads\t" total/n;
         print "Avg_Mapped_Reads\t" mapped/n;
         print "Avg_Percent_Mapped\t" pct/n;
     }' $OUT >> $OUT

echo "Done."
