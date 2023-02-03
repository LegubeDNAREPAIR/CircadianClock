## bigwig file path with inverted commas
## you can do more than one file at once on separate plots by doing: "file1" "file"
## the same for the sample and plot labels: "lab1" "lab2"

## If you just want to have more than one bigwigs on the same plot then you can
## put the files in one set of inverted commas separated by a space e.g. "file1 file2"
## You will then need to do the same for the sample labels (LABS) but not for the plot "lab1 lab2"
## labels (NAMES), as there is still only one name per plot


declare -a FILES=(
"bigwigCompare_HC2FWBGXH_PER2OHTvsPER2DIVA_log2_bs50.bw"
"bigwigCompare_HC2FWBGXH_BMAL1OHTvsBMAL1DIVA_log2_bs50.bw"
)


## sample label, one per bigwig
declare -a LABS=(
"bwC_PER2"
"bwC_BMAL1"
)

## plot label, one per plot 
declare -a NAMES=(
"bwC_PER2"
"bwC_BMAL1"
)

## the same rules apply if you want to bed files to be on the same plot

## bed file path
declare -a BED=(
"BLESS_HR_JunFragPE_Rmdups_pm500bp.bed"
"BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed"
"BLESS_80best_JunFragPE_Rmdups_pm500bp.bed"
"80random.bed"
)


## bed file label
declare -a BED_LABS=(
"HR_DSB"
"NHEJ_DSB"
"80_DSB"
"80_random"
)


## file regions name (to be used in the plot title)
declare -a FILE_EXT=(
"HR_DSB" 
"NHEJ_DSB"
"80_DSB" 
"80_random"
)

## total region to plot from center of file
## write the total region length for plot title, and then the half length in bp for the command itself
KB="10kb"
KB_bp=5000

## bin size - i usually try to keep the profiles to have 200-400 bins
bs=50

## out file path
OUT="/media/scollins/Sarah_scRNA/PER2_BMAL1/"


###########################################################
# you dont need to change anything after this line
###########################################################


nrow_bed=$((${#BED[@]}-1))
nrow=$((${#NAMES[@]}-1))
echo $nrow


for bed in $(seq 0 $nrow_bed)
do
FILE_EXT=$(echo "${FILE_EXT[$bed]}")
BED_LABS=$(echo "${BED_LABS[$bed]}")
BED=$(echo "${BED[$bed]}")

echo $FILE_EXT
echo $BED_LABS
echo $BED

for i in $(seq 0 $nrow)
do
FILE=$(echo "${FILES[$i]}")
LAB=$(echo "${LABS[$i]}")
NAME=$(echo "${NAMES[$i]}")

echo $FILE
echo $LAB
echo $NAME

echo " "
echo "Running...."

## -a and -b refer to base pairs to plot from the center
computeMatrix reference-point -S ${FILE} \
-R ${BED} \
--referencePoint center -a $KB_bp -b $KB_bp \
--samplesLabel $LAB \
--outFileName ${OUT}${NAME}_${FILE_EXT}_${KB}.ma.gz \
--binSize $bs \
--sortRegions keep \
-p 15

echo "Done!"
echo " "

## Only use this if plotting over genes bodies
#computeMatrix scale-regions -S ${FILE} \
#-R ${BED} \
#-a 2000 -b 2000 \
#--samplesLabel $LAB \
#--regionBodyLength 5000 \
#--outFileName ${OUT}${NAME}_${FILE_EXT}_${KB}.ma.gz \
#--binSize 50 \
#--sortRegions descend \
#-p 20

done
done
wait





