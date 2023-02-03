## APA plot to process output of APA from juicertools


require(tidyverse)
require(plyranges)
require(reshape2)
require(rtracklayer)

seqlens <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19 %>% seqlengths()
#Functions

### process APA plots

process_APA <- function(file){
  seqpos <- seq(-10,10,by=1) %>% as.character()
  my.dat <- file %>% read_csv(col_names = F) %>%
    mutate_if(is.character,str_remove,"\\[|\\]") %>%
    mutate_if(is.character,as.numeric) %>%
    as.matrix()
  
  colnames(my.dat) <- seqpos
  rownames(my.dat) <- rev( seqpos )
  my.dat <- my.dat %>% reshape2::melt() %>%
    mutate(file = file) 
  return(my.dat)
}


files <- c(
  "siCTRL_DIVA_80DSB"="/home/scollins/Hic/siCTRL_DIVA_BEN/Intrachrom_noTads_BLESS_80_nochr.bed/10000/gw/APA.txt",
  "siCTRL_OHT_80DSB"="/home/scollins/Hic/siCTRL_OHT_BEN/Intrachrom_noTads_BLESS_80_nochr.bed/10000/gw/APA.txt",
  "siPER2_DIVA_80DSB"="/home/scollins/Hic/siPER2_DIVA_BEN/Intrachrom_noTads_BLESS_80_nochr.bed/10000/gw/APA.txt",
  "siPER2_OHT_80DSB"="/home/scollins/Hic/siPER2_OHT_BEN/Intrachrom_noTads_BLESS_80_nochr.bed/10000/gw/APA.txt",
  "siCTRL_DIVA_80random"="/home/scollins/Hic/siCTRL_DIVA_BEN/Intrachrom_withTads_random_80_nochr.bed/10000/gw/APA.txt",
  "siCTRL_OHT_80random"="/home/scollins/Hic/siCTRL_OHT_BEN/Intrachrom_withTads_random_80_nochr.bed/10000/gw/APA.txt",
  "siPER2_DIVA_80random"="/home/scollins/Hic/siPER2_DIVA_BEN/Intrachrom_withTads_random_80_nochr.bed/10000/gw/APA.txt",
  "siPER2_OHT_80random"="/home/scollins/Hic/siPER2_OHT_BEN/Intrachrom_withTads_random_80_nochr.bed/10000/gw/APA.txt"
) 




dat.apa <- NULL
dat.fc <- NULL
for(i in 1:length(files)){
  df <- process_APA(files[i])
  df$Condition <- names(files)[i]
  
  mid_val <- subset(df, Var1 %in% c(0,1,-1) & Var2 %in% c(0,1,-1))
  mid_val <- mean(mid_val$value)
  
  low_val <- subset(df, Var1 %in% c(-10,-9,-8) & Var2 %in% c(-10,-9,-8))
  low_val <- mean(low_val$value)
  
  fc = mid_val / low_val
  
  df.fc <- data.frame(Condition=names(files)[i],
                      FC=fc)
  
  dat.apa <- rbind(dat.apa,df)
  dat.fc <- rbind(dat.fc,df.fc)
}

dat.apa$Condition <- factor(dat.apa$Condition, levels = names(files))

dat.apa <- dat.apa %>% left_join(dat.fc)
dat.apa$FC <- format(round(dat.apa$FC, 2), nsmall = 2)
dat.apa$label <- paste0(dat.apa$Condition," \n FC = ",dat.apa$FC)

res="100kb"
res1="10kb"
title="HiC siCTRL siPER2 80DSB and CTRL -/+ 100kb,"


p1 <- dat.apa %>% ggplot(aes(Var2,Var1)) + 
  geom_raster(aes(fill = value, color=value)) +
  scale_fill_gradient2(low = "white",high = "#EA2027"
  )  +
  scale_color_gradient2(low = "white",high = "#EA2027") +
  theme_classic(base_size = 8) +
  facet_wrap(~label,nrow=3, scales = "free") +
  ylab("") + xlab("") +
  scale_x_continuous(breaks=c(-10,0,10),
                     labels=c(paste0("-",res), "DSB", paste0("+",res))) +
  scale_y_continuous(breaks=c(-10,0,10),
                     labels=c(paste0("-",res), "DSB", paste0("+",res))) +
  ggtitle(paste0(title," resolution: ", res1))

p1
