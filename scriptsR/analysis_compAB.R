# Raphael Mourad
# 07/02/2019


# GOOD RESULTS!



# WORKING DIRECTORY ----------------------------------------------------

setwd("/media/mourad/diskSave/MCF_Toulouse/recherche/LegubeTeam/")

# LIBRARIES ------------------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg19)
require(tidyverse)
require(plyranges)


# Genome 
seqInfohg19=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
seqlen=seqlengths(seqInfohg19)

# Parameters
binSize=100000
binSizeTxt="100kb"

# Chromosomes
SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
Chr.V=glue::glue("chr{c(1:22,'X')}")
seqlen=seqlengths(SeqInfo[Chr.V])

options(scipen=999)

# ATACseq
ATAC.GR=read_narrowpeaks("/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ATAC-Seq/Clouaire_HC3HCBGX9-A-BColine/PROCESSED/mapping/mapQ25/bam/HC3HCBGX9_ATACseqA_DIvA_18s005247-1-1_Clouaire_lane1ATACseqADIvA/HC3HCBGX9_ATACseqA_DIvA_18s005247-1-1_Clouaire_lane1ATACseqADIvA_peaks.narrowPeak")


files.comp <- list.files("AB/KR",full.names = T)
# combi_files <- Chr.V %>% map(function(x){str_subset(files.comp,pattern=glue::glue("_{x}_.+_{binSizeTxt}"))})
combi_files <- setNames(c("siCTRL_DIVA","siCTRL_OHT","siPER2_DIVA","siPER2_OHT"),c("siCTRL_DIVA","siCTRL_OHT","siPER2_DIVA","siPER2_OHT")) %>%
  map(function(x){str_subset(files.comp,pattern=glue::glue("_{x}_{binSizeTxt}"))})

# Load compAB
for(i in names(combi_files)){
 
 filesAB <-combi_files[[i]]
 bini.GR= tileGenome(seqlen[Chr.V],tilewidth = binSize,cut.last.tile.in.chrom = T)
 res <- filesAB %>% lapply(function(x){
   # message(x)
   score <- read_tsv(x,col_names = F) %>% pull(X1) %>% scale(scale=T) %>% as.vector()
   sres <- bini.GR %>% filter(seqnames == str_extract(x,"chr[0-9X]+")) %>% mutate(score = score)
   ATACbini=countOverlaps(sres,ATAC.GR)
   if(cor(sres$score,ATACbini,use="pairwise.complete.obs")<0){
     sres$score <- -1*sres$score
   }
   sres$score[is.na(sres$score)] <- 0
   return(sres)
 })
 res <- do.call("c",res)
 write_bed_graph(res,glue::glue("AB/compartimentAB_KR_{i}_{binSizeTxt}.bedGraph"))

}
# - mettre les data sur le NAS + lien cloud 
# - verifier les compD sur le NAS si ils y sont 
# figurez S2D du papier -> normalement s'est fait mais verifier la windows
# - jouer sur fenêtre / contraste des boxplots
# - tester pm 40 kb
# - faire compartiments A et B de la manip de Benj
# - saddle plot Per2 sur Control
# - S1D -> within damaged TAD contact DSB + Random


