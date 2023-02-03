
## APA plots from juicertools

JUICER="java -jar /home/scollins/Tools/Juicer/scripts/juicer_tools_1.22.01.jar"

## HiC files
declare -a FILES=(
"/media/nas1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Hi-C/HiC_F_Ben/siCTRL_DIVA/aligned/inter.hic"
"/media/nas1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Hi-C/HiC_F_Ben/siCTRL_OHT/aligned/inter.hic"
"/media/nas1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Hi-C/HiC_F_Ben/siPER2_DIVA/aligned/inter.hic"
"/media/nas1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Hi-C/HiC_F_Ben/siPER2_OHT/aligned/inter.hic"
)

## HiC project name
declare -a LABS=(
"siCTRL_DIVA_BEN"
"siCTRL_OHT_BEN"
"siPER2_DIVA_BEN"
"siPER2_OHT_BEN"
)

## Bed files
declare -a BED=(
"/home/scollins/Hic/Interchrom_80_random_nochr.bed"
"/home/scollins/Hic/Interchrom_80_DSB_nochr.bed"
"/home/scollins/Hic/Intrachrom_80_random_nochr.bed"
"/home/scollins/Hic/Intrachrom_80_DSB_nochr.bed")

## Bed file labels
declare -a BED_LABS=(
"Interchrom_80_random"
"Interchrom_80_DSB"
"Intrachrom_80_random"
"Intrachrom_80_DSB"
)

nrow_bed=$((${#BED[@]}-1))
nrow=$((${#LABS[@]}-1))
echo $nrow


for bed in $(seq 0 $nrow_bed)
do
BED_LABS=$(echo "${BED_LABS[$bed]}")
BED=$(echo "${BED[$bed]}")

echo $BED_LABS
echo $BED

for i in $(seq 0 $nrow)
do
FILE=$(echo "${FILES[$i]}")
LAB=$(echo "${LABS[$i]}")

mkdir $LAB
echo $FILE
echo $LAB

$JUICER apa -r 10000 $FILE $BED ${LAB}/${BED_LABS} --threads 30 -e

done
done
