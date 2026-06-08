#!/bin/bash
set -euo pipefail

BASE="metadata/superpops"

# Create superpopulation directories
mkdir -p $BASE/AFR $BASE/EUR $BASE/EAS $BASE/SAS

# Helper function
add_entry () {
    local superpop=$1
    local pop=$2
    local replacement=$3
    local url1=$4
    local url2=$5

    echo -e "${replacement}\t${url1}\t${url2}" >> "${BASE}/${superpop}/${pop}.txt"
}

#############################################
# AFR — ACB
#############################################
add_entry AFR ACB HG01880 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047764/ERR047764_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047764/ERR047764_2.fastq.gz
add_entry AFR ACB HG01882 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047765/ERR047765_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047765/ERR047765_2.fastq.gz
add_entry AFR ACB HG01886 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR360/SRR360655/SRR360655_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR360/SRR360655/SRR360655_2.fastq.gz
add_entry AFR ACB HG02013 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047788/ERR047788_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047788/ERR047788_2.fastq.gz
add_entry AFR ACB HG01890 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR710/SRR710118/SRR710118_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR710/SRR710118/SRR710118_2.fastq.gz

#############################################
# AFR — ASW
#############################################
add_entry AFR ASW NA20274 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250968/ERR250968_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250968/ERR250968_2.fastq.gz
add_entry AFR ASW NA19922 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR708/SRR708373/SRR708373_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR708/SRR708373/SRR708373_2.fastq.gz
add_entry AFR ASW NA20342 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR071/SRR071187/SRR071187_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR071/SRR071187/SRR071187_2.fastq.gz
add_entry AFR ASW NA20351 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR708/SRR708375/SRR708375_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR708/SRR708375/SRR708375_2.fastq.gz
add_entry AFR ASW NA20355 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250976/ERR250976_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250976/ERR250976_2.fastq.gz

#############################################
# AFR — ESN
#############################################
add_entry AFR ESN HG02977 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250481/ERR250481_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250481/ERR250481_2.fastq.gz
add_entry AFR ESN HG02981 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250485/ERR250485_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250485/ERR250485_2.fastq.gz
add_entry AFR ESN HG03157 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250541/ERR250541_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250541/ERR250541_2.fastq.gz
add_entry AFR ESN HG03139 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250539/ERR250539_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250539/ERR250539_2.fastq.gz
add_entry AFR ESN HG03267 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250577/ERR250577_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250577/ERR250577_2.fastq.gz

#############################################
# AFR — LWK
#############################################
add_entry AFR LWK NA19020 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR748/SRR748550/SRR748550_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR748/SRR748550/SRR748550_2.fastq.gz
add_entry AFR LWK NA19028 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR748/SRR748773/SRR748773_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR748/SRR748773/SRR748773_2.fastq.gz

#############################################
# AFR — YRI
#############################################
add_entry AFR YRI NA18488 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100014/SRR100014_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100014/SRR100014_2.fastq.gz
add_entry AFR YRI NA18498 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100024/SRR100024_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100024/SRR100024_2.fastq.gz
add_entry AFR YRI NA18499 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100023/SRR100023_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100023/SRR100023_2.fastq.gz
add_entry AFR YRI NA18870 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100031/SRR100031_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR100/SRR100031/SRR100031_2.fastq.gz

#############################################
# EUR — FIN
#############################################
add_entry EUR FIN HG00171 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034564/ERR034564_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034564/ERR034564_2.fastq.gz
add_entry EUR FIN HG00173 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/SRR764756/SRR764756_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR764/SRR764756/SRR764756_2.fastq.gz

#############################################
# EUR — IBS
#############################################
add_entry EUR IBS HG01501 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR701/SRR701470/SRR701470_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR701/SRR701470/SRR701470_2.fastq.gz
add_entry EUR IBS HG01503 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR709/SRR709962/SRR709962_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR709/SRR709962/SRR709962_2.fastq.gz
add_entry EUR IBS HG01504 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR710/SRR710112/SRR710112_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR710/SRR710112/SRR710112_2.fastq.gz

#############################################
# EUR — GBR
#############################################
add_entry EUR GBR HG00097 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR765/SRR765989/SRR765989_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR765/SRR765989/SRR765989_2.fastq.gz
add_entry EUR GBR HG00105 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR711/SRR711354/SRR711354_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR711/SRR711354/SRR711354_2.fastq.gz
add_entry EUR GBR HG00106 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR702/SRR702073/SRR702073_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR702/SRR702073/SRR702073_2.fastq.gz
add_entry EUR GBR HG00125 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR031/ERR031932/ERR031932_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR031/ERR031932/ERR031932_2.fastq.gz
add_entry EUR GBR HG00145 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099957/SRR099957_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099957/SRR099957_2.fastq.gz

#############################################
# EUR — TSI
#############################################
add_entry EUR TSI NA20503 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR766/SRR766036/SRR766036_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR766/SRR766036/SRR766036_2.fastq.gz
add_entry EUR TSI NA20511 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250978/ERR250978_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250978/ERR250978_2.fastq.gz
add_entry EUR TSI NA20507 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR748/SRR748292/SRR748292_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR748/SRR748292/SRR748292_2.fastq.gz

#############################################
# EAS — CDX
#############################################
add_entry EAS CDX HG00956 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR050/ERR050740/ERR050740_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR050/ERR050740/ERR050740_2.fastq.gz

#############################################
# EAS — CHB
#############################################
add_entry EAS CHB NA18530 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034576/ERR034576_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR034/ERR034576/ERR034576_2.fastq.gz

#############################################
# EAS — CHS
#############################################
add_entry EAS CHS HG00407 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR031/ERR031987/ERR031987_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR031/ERR031987/ERR031987_2.fastq.gz

#############################################
# EAS — JPT
#############################################
add_entry EAS JPT NA18946 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR718/SRR718079/SRR718079_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR718/SRR718079/SRR718079_2.fastq.gz
add_entry EAS JPT NA19001 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR707/SRR707200/SRR707200_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR707/SRR707200/SRR707200_2.fastq.gz

#############################################
# EAS — KHV
#############################################
add_entry EAS KHV HG01596 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047699/ERR047699_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047699/ERR047699_2.fastq.gz

#############################################
# SAS — GIH
#############################################
add_entry SAS GIH NA20864 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047869/ERR047869_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR047/ERR047869/ERR047869_2.fastq.gz
add_entry SAS GIH NA21092 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR360/SRR360498/SRR360498_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR360/SRR360498/SRR360498_2.fastq.gz

#############################################
# SAS — PJL
#############################################
add_entry SAS PJL HG02651 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250389/ERR250389_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250389/ERR250389_2.fastq.gz
add_entry SAS PJL HG02691 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250399/ERR250399_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250399/ERR250399_2.fastq.gz

#############################################
# SAS — ITU
#############################################
add_entry SAS ITU HG03723 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250731/ERR250731_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250731/ERR250731_2.fastq.gz
add_entry SAS ITU HG03790 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250761/ERR250761_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250761/ERR250761_2.fastq.gz

#############################################
# SAS — STU
#############################################
add_entry SAS STU HG03646 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250689/ERR250689_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250690/ERR250690_1.fastq.gz
add_entry SAS STU HG03672 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250695/ERR250695_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250695/ERR250695_2.fastq.gz

#############################################
# SAS — BEB
#############################################
add_entry SAS BEB HG03603 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250671/ERR250671_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250671/ERR250671_2.fastq.gz
add_entry SAS BEB HG03814 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250773/ERR250773_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250773/ERR250773_2.fastq.gz
add_entry SAS BEB HG04144 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250872/ERR250872_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR250/ERR250872/ERR250872_2.fastq.gz
