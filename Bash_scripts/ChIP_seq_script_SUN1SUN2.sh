GENOME="/bioinfo_legube/VINCENT/GENOMES/female.hg19/female.hg19.fa"
DIR="/bioinfo_legube/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Clouaire_AACJ327M5_CTCF_24H_SUN1SUN2/"
PREFIX="AACJ327M5_Pool_100123_23s000062-1-1_Clouaire_lane1"

declare -a NAMES=(
"SUN1DIvA"
"SUN1OHT"
"SUN2DIvA"
"SUN2OHT"
)

nrow=$((${#NAMES[@]} -1))

mkdir ${DIR}mapping
mkdir ${DIR}mapping/bam/
mkdir ${DIR}mapping/stats/
mkdir ${DIR}mapping/bigwig/

for i in $(seq 0 $nrow)
do
NAME=$(echo "${NAMES[$i]}")

## align and direct output to sort
bwa mem -t 30 $GENOME $DIR/RAW/${PREFIX}${NAME}_sequence.txt.gz | samtools sort -@ 24 -o $DIR/mapping/bam/${PREFIX}${NAME}_sorted.bam -

## rmdups
samtools rmdup -s $DIR/mapping/bam/${PREFIX}${NAME}_sorted.bam $DIR/mapping/bam/${PREFIX}${NAME}_rmdups.bam

## Samtools index
samtools index  $DIR/mapping/bam/${PREFIX}${NAME}_rmdups.bam $DIR/mapping/bam/${PREFIX}${NAME}_rmdups.bai

samtools stats $DIR/mapping/bam/${PREFIX}${NAME}_rmdups.bam > $DIR/mapping/stats/${PREFIX}${NAME}_rmdups.stats

## Samtools idxstats
samtools idxstats $DIR/mapping/bam/${PREFIX}${NAME}_rmdups.bam >  $DIR/mapping/stats/${PREFIX}${NAME}_rmdups.idxstats

## Samtools flagstat
samtools flagstat  $DIR/mapping/bam/${PREFIX}${NAME}_rmdups.bam >  $DIR/mapping/stats/${PREFIX}${NAME}_rmdups.flagstat

bamCoverage -b $DIR/mapping/bam/${PREFIX}${NAME}_rmdups.bam -o $DIR/mapping/bigwig/${PREFIX}${NAME}_normalized_1.bw --exactScaling --normalizeUsing CPM -bs 1 -p 30

done
