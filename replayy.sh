#!/bin/bash
#SBATCH --job-name=dl_missing_repl
#SBATCH --output=logs/dl_missing_%A_%a.out
#SBATCH --error=logs/dl_missing_%A_%a.err
#SBATCH --array=0-13%14
#SBATCH --cpus-per-task=5
#SBATCH --mem=14G

set -euo pipefail
mkdir -p logs

#############################################
# HARD‑CODED ARRAYS FOR MISSING REPLACEMENTS
#############################################

SUPERPOP=(
AFR
AFR
AFR
AFR
AFR
AFR
AFR
EUR
EUR
EUR
EAS
EAS
EAS
SAS
)

POP=(
GWD
GWD
GWD
GWD
GWD
YRI
YRI
CEU
TSI
TSI
JPT
JPT
KHV
BEB
)


SAMPLEID=(
HG02813
HG02814
HG02840
HG02878
HG02982
NA18499
NA18870
NA06989
NA20511
NA20507
NA18946
NA19001
HG01596
HG04144
)

URL1=(
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250421/ERR250421_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250423/ERR250423_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250435/ERR250435_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250445/ERR250445_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250487/ERR250487_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100023/SRR100023_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100031/SRR100031_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/SRR098417/SRR098417_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250978/ERR250978_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR748/SRR748292/SRR748292_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR718/SRR718079/SRR718079_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR707/SRR707200/SRR707200_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047699/ERR047699_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250872/ERR250872_1.fastq.gz
)

URL2=(
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250421/ERR250421_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250423/ERR250423_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250435/ERR250435_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250445/ERR250445_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250487/ERR250487_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100023/SRR100023_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100031/SRR100031_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR098/SRR098417/SRR098417_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250978/ERR250978_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR748/SRR748292/SRR748292_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR718/SRR718079/SRR718079_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR707/SRR707200/SRR707200_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047699/ERR047699_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250872/ERR250872_2.fastq.gz
)

#############################################
# SELECT SAMPLE FOR THIS ARRAY TASK
#############################################

IDX=$SLURM_ARRAY_TASK_ID

SP=${SUPERPOP[$IDX]}
POP=${POP[$IDX]}
SID=${SAMPLEID[$IDX]}
R1=${URL1[$IDX]}
R2=${URL2[$IDX]}

OUTDIR="data/${SP}/${POP}/${SID}"
mkdir -p "$OUTDIR"

echo "[$(date)] Downloading $SID into $OUTDIR"

wget -c "$R1" -O "${OUTDIR}/${SID}_1.fastq.gz"
wget -c "$R2" -O "${OUTDIR}/${SID}_2.fastq.gz"

echo "[$(date)] Finished $SID"
