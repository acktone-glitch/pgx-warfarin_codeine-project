#!/bin/bash
#SBATCH --job-name=dl_replacements
#SBATCH --output=logs/dl_replacements_%A_%a.out
#SBATCH --error=logs/dl_replacements_%A_%a.err
#SBATCH --array=0-50%10
#SBATCH --cpus-per-task=64
#SBATCH --mem=48G

set -euo pipefail

mkdir -p logs

# One entry per replacement sample, same index across all arrays
SUPERPOP=(
AFR AFR AFR AFR AFR
AFR AFR AFR AFR AFR
AFR AFR AFR AFR AFR
AFR AFR AFR AFR
EUR EUR
EUR EUR EUR
EUR EUR EUR EUR EUR
EUR EUR EUR
EAS EAS EAS EAS EAS
SAS SAS
SAS SAS
SAS SAS
SAS SAS
SAS SAS
SAS SAS SAS
)

POP=(
ACB ACB ACB ACB ACB
ASW ASW ASW ASW ASW
ESN ESN ESN ESN ESN
LWK LWK
YRI YRI YRI YRI
FIN FIN
IBS IBS IBS
GBR GBR GBR GBR GBR
TSI TSI TSI
CDX CHB CHS JPT JPT
KHV
GIH GIH
PJL PJL
ITU ITU
STU STU
BEB BEB BEB
)

SAMPLEID=(
HG01880
HG01882
HG01886
HG02013
HG01890
NA20274
NA19922
NA20342
NA20351
NA20355
HG02977
HG02981
HG03157
HG03139
HG03267
NA19020
NA19028
NA18488
NA18498
NA18499
NA18870
HG00171
HG00173
HG01501
HG01503
HG01504
HG00097
HG00105
HG00106
HG00125
HG00145
NA20503
NA20511
NA20507
HG00956
NA18530
HG00407
NA18946
NA19001
HG01596
NA20864
NA21092
HG02651
HG02691
HG03723
HG03790
HG03646
HG03672
HG03603
HG03814
HG04144
)

URL1=(
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047764/ERR047764_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047765/ERR047765_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR360/SRR360655/SRR360655_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047788/ERR047788_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR710/SRR710118/SRR710118_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250968/ERR250968_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR708/SRR708373/SRR708373_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR071/SRR071187/SRR071187_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR708/SRR708375/SRR708375_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250976/ERR250976_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250481/ERR250481_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250485/ERR250485_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250541/ERR250541_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250539/ERR250539_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250577/ERR250577_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR748/SRR748550/SRR748550_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR748/SRR748773/SRR748773_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100014/SRR100014_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100024/SRR100024_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100023/SRR100023_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100031/SRR100031_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034564/ERR034564_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/SRR764756/SRR764756_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR701/SRR701470/SRR701470_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR709/SRR709962/SRR709962_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR710/SRR710112/SRR710112_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR765/SRR765989/SRR765989_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR711/SRR711354/SRR711354_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR702/SRR702073/SRR702073_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR031/ERR031932/ERR031932_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099957/SRR099957_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR766/SRR766036/SRR766036_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250978/ERR250978_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR748/SRR748292/SRR748292_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR050/ERR050740/ERR050740_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034576/ERR034576_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR031/ERR031987/ERR031987_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR718/SRR718079/SRR718079_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR707/SRR707200/SRR707200_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047699/ERR047699_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047869/ERR047869_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR360/SRR360498/SRR360498_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250389/ERR250389_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250399/ERR250399_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250731/ERR250731_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250761/ERR250761_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250689/ERR250689_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250695/ERR250695_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250671/ERR250671_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250773/ERR250773_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250872/ERR250872_1.fastq.gz
)

URL2=(
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047764/ERR047764_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047765/ERR047765_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR360/SRR360655/SRR360655_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047788/ERR047788_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR710/SRR710118/SRR710118_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250968/ERR250968_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR708/SRR708373/SRR708373_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR071/SRR071187/SRR071187_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR708/SRR708375/SRR708375_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250976/ERR250976_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250481/ERR250481_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250485/ERR250485_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250541/ERR250541_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250539/ERR250539_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250577/ERR250577_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR748/SRR748550/SRR748550_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR748/SRR748773/SRR748773_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100014/SRR100014_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100024/SRR100024_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100023/SRR100023_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100031/SRR100031_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034564/ERR034564_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/SRR764756/SRR764756_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR701/SRR701470/SRR701470_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR709/SRR709962/SRR709962_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR710/SRR710112/SRR710112_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR765/SRR765989/SRR765989_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR711/SRR711354/SRR711354_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR702/SRR702073/SRR702073_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR031/ERR031932/ERR031932_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099957/SRR099957_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR766/SRR766036/SRR766036_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250978/ERR250978_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR748/SRR748292/SRR748292_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR050/ERR050740/ERR050740_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034576/ERR034576_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR031/ERR031987/ERR031987_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR718/SRR718079/SRR718079_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR707/SRR707200/SRR707200_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047699/ERR047699_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047869/ERR047869_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR360/SRR360498/SRR360498_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250389/ERR250389_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250399/ERR250399_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250731/ERR250731_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250761/ERR250761_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250690/ERR250690_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250695/ERR250695_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250671/ERR250671_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250773/ERR250773_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250872/ERR250872_2.fastq.gz
)

IDX=${SLURM_ARRAY_TASK_ID}

SP=${SUPERPOP[$IDX]}
POP=${POP[$IDX]}
SID=${SAMPLEID[$IDX]}
R1=${URL1[$IDX]}
R2=${URL2[$IDX]}

OUTDIR="data/${SP}/${POP}/${SID}"
mkdir -p "$OUTDIR"

echo "[$(date)] Downloading ${SID} into ${OUTDIR}"

wget -c "$R1" -O "${OUTDIR}/${SID}_1.fastq.gz"
wget -c "$R2" -O "${OUTDIR}/${SID}_2.fastq.gz"

echo "[$(date)] Done ${SID}"
