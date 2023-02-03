# require(HiTC)
# require(tidyverse)
# require(plyranges)
# require(data.table)
# require(Matrix)
# require(BSgenome.Hsapiens.UCSC.hg19)
# DSB174 <- read_bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/ASIsites_hg19_174_clived_IQ1.5.bed")
# loadRData <- function(fileName){
#   #loads an RData file, and returns it
#   load(fileName)
#   get(ls()[ls() != "fileName"])
# }
# 
# read_Dumped_matrix <- function(fileCounti,binendi = NULL){
#   datai=as.matrix(fread(paste0("zcat ",fileCounti),sep='\t',header=F))
#   
#   datai=datai[!is.na(datai[,3]),]
#   datai=rbind(datai,c((binendi-1)*binSize,(binendi-1)*binSize,0))
#   if(obserOE=="observed"){
#     datai[,3]=round(datai[,3])
#   }else{
#     datai[,3]=round(datai[,3],2)
#   }
#   return(datai)
# }
# 
# obserOE <- "OE"
# # expe <- c("HiC_D_DIvA","HiC_D_OHT","HiC_D_OHTATMi","HiC_D_OHTDNAPKi","HiC_D_OHTPARPi")
# expe <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_BENJ/data/",pattern=".hic") %>% str_remove(".hic")
# # dir_dump <- "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/"
# dir_dump <- "/home/rochevin/Documents/PROJET_INGE/HiC_Aline/dump/TRANS/KR"
# # Resolution
# binSize=100000
# binSizeTxt=paste0(binSize/1e3,"kb")
# 
# SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
# 
# my_chromosomes <- seqlevels(DSB174)
# 
# DIvA <- expe[1]
# # OHT <- expe[4]
# for(OHT in expe[c(2)]){
#   res_pca <- lapply(my_chromosomes,function(Chr.V){
#     message(Chr.V)
#     chrendi=seqlengths(SeqInfo[Chr.V])
#     binendi=ceiling(chrendi/binSize)
#     
#     fileCounti_DIvA=glue::glue("{dir_dump}/dump_{obserOE}_KR_{DIvA}_{Chr.V}_{Chr.V}_{binSizeTxt}.txt.gz")
#     # fileCounti_OHT=paste0("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/data/KR/dump_",obserOE,"_KR_",expe[3],"_",Chr.V,"_",Chr.V,"_",binSizeTxt,".txt.gz")
#     fileCounti_OHT=glue::glue("{dir_dump}/dump_{obserOE}_KR_{OHT}_{Chr.V}_{Chr.V}_{binSizeTxt}.txt.gz")
#     
#     
#     
#     data_mOHT <- read_Dumped_matrix(fileCounti_DIvA,binendi=binendi) %>% 
#       as_tibble %>% 
#       setNames(c("bin1","bin2","mOHT"))
#     data_pOHT <- read_Dumped_matrix(fileCounti_OHT,binendi=binendi) %>% 
#       as_tibble() %>%
#       setNames(c("bin1","bin2","pOHT")) %>% inner_join(data_mOHT,by = c("bin1","bin2"))
#     
#     datai <- data_pOHT %>% mutate(ratio = log2(pOHT/mOHT))
#     datai <- datai %>% mutate(ratio = ifelse(pOHT==0 & mOHT== 0,0,ratio))
#     datai <- datai[is.finite(datai$ratio),]
#     data.Mati=sparseMatrix(i=(datai$bin1/binSize)+1,j=(datai$bin2/binSize)+1,x=datai$ratio)
#     data.MatSymi=data.Mati+t(data.Mati)
#     diag(data.MatSymi)=diag(data.MatSymi)/2
#     rm(data.Mati,datai)
#     
#     starti=seq(1,chrendi,by=binSize)
#     endi=c(seq(binSize,chrendi,by=binSize),chrendi)
#     chri.GR=GRanges(Chr.V,IRanges(starti,endi))
#     names(chri.GR)=paste0("bin",1:length(chri.GR))
#     
#     HTC=HTCexp(data.MatSymi,chri.GR,chri.GR)
#     
#     cm <- HiTC:::sparseCor(intdata(HTC))
#     cm[is.na(cm)] <- 0
#     princp = princomp(cm)
#     
#     
#     PCA_1 <-bind_cols(as_tibble(x_intervals(HTC)),as_tibble(loadings(princp)[])) %>%
#       as_granges()%>% join_overlap_left(DSB174) 
#     PCA_cor_DSB <- PCA_1 %>% mutate(DSB = ifelse(!is.na(name),1,0)) %>% as_tibble() %>% dplyr::select(-seqnames:-strand,-name,-score)
#     
#     
#     PCA_cor_DSB <- PCA_cor_DSB %>% dplyr::select(-DSB) %>% apply(2,function(one_col){
#       cc <- cor.test(one_col,PCA_cor_DSB$DSB)
#       tibble(my_cor= cc$estimate,p.val = cc$p.value)
#     }) %>% bind_rows() %>% mutate(name =colnames(PCA_cor_DSB)[-length(PCA_cor_DSB)])
#     
#     
#     
#     PCA_1 <- PCA_1 %>% as_tibble() %>%
#       gather(key = PCs,value = value,-seqnames:-strand,-name:-score) %>% 
#       left_join(PCA_cor_DSB,by = c("PCs"="name"))
#     
#     PCA_1 %>% filter(PCs %in% str_c("Comp",1:10,sep = ".")) 
#     
#   }) %>% bind_rows()
#   
#   
#   
#   
#   res_pca <- res_pca  %>% as_granges()
#   
#   # cc1 <- res_pca %>% filter_by_overlaps(DSB174 %>% anchor_center() %>% mutate(width = 10000)) %>% as_tibble() %>% group_by(seqnames) %>%
#   #   summarise(meanval = mean(value)) %>% filter(meanval < 0) %>% pull(seqnames) %>% as.character()
#   # 
#   # 
#   # res_pca <- res_pca %>% as_tibble() %>%split(.,.$seqnames)
#   # res_pca <- res_pca[-1]
#   # res_pca <- res_pca %>% map(function(x){
#   #   if(unique(x$seqnames) %in% cc1){
#   #     x %>% mutate(value = value * -1)
#   #   }else{
#   #     x
#   #   }
#   # })  
#   
#   res_pca <- res_pca  %>% as_granges()
#   res_pca <- res_pca %>% as_tibble() %>%split(.,.$seqnames)
#   res_pca <- res_pca[-1]
#   
#   res_p <- lapply(res_pca[my_chromosomes],function(one_data){
#     one_data <- one_data %>% filter(PCs %in% str_c("Comp",1,sep = "."))%>% mutate(my_cor = round(my_cor,3))%>% mutate(p.val = format.pval(p.val,3))%>% unite("PCs",PCs,my_cor,sep = " cor: ") %>% unite("PCs",PCs,p.val,sep = " p.val: ") 
#     one_data %>%  ggplot(aes(x=start,y=value)) + geom_line() +
#       facet_wrap(~PCs,scales="free",ncol=1) +
#       geom_vline(data = filter(one_data,!is.na(name)),aes(xintercept=start),linetype ="dashed") +ggtitle(unique(one_data$seqnames))
#   })
#   
#   pdf(glue::glue("results/{OHT}_test_pca_cor_alamano_100kb_1_with_reverse_signal_{binSizeTxt}.pdf"),height=4,width=12)
#   print(res_p)
#   dev.off()
#   
#   ##Extract to bigwig 
#   cc <- res_pca %>% map(filter,PCs %in% str_c("Comp",1,sep = ".")) %>% map(mutate,value=ifelse(my_cor < 0,value*-1,value)) %>% 
#     bind_rows() %>% dplyr::select(seqnames:strand,value) %>% as_granges()
#   # cc <- res_pca %>% bind_rows() %>% filter(PCs == "Comp.1") %>% dplyr::select(seqnames:strand,value) %>% as_granges()
#   
#   
#   setwd("/media/HDD_ROCHER/PROJET_INGE/HiC_BENJ/")
#   cc <- regioneR::filterChromosomes(cc,keep.chr=my_chromosomes)
#   seqlengths(cc) <- seqlengths(SeqInfo)[names(seqlengths(cc))]
#   cc.cov <- coverage(cc,weight="value")
#   export.bw(cc.cov,glue::glue("results/PC1_all_chr_log2ratio_{binSizeTxt}_{OHT}.bw"))
#   
# }

##V2 WITH WHAT GAELLE ASKED

require(HiTC)
require(tidyverse)
require(plyranges)
require(data.table)
require(Matrix)
require(BSgenome.Hsapiens.UCSC.hg19)
DSB174 <- read_bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/ASIsites_hg19_174_clived_IQ1.5.bed")
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

read_Dumped_matrix <- function(fileCounti,binendi = NULL){
  datai=as.matrix(fread(paste0("zcat ",fileCounti),sep='\t',header=F))
  
  datai=datai[!is.na(datai[,3]),]
  datai=rbind(datai,c((binendi-1)*binSize,(binendi-1)*binSize,0))
  if(obserOE=="observed"){
    datai[,3]=round(datai[,3])
  }else{
    datai[,3]=round(datai[,3],2)
  }
  return(datai)
}

Get1val <- function(my.wigs,one.w,x){
  lapply(split(x,droplevels(seqnames(x))),function(zz){
    message(unique(as.character(seqnames(zz))))
    cov <- one.w[[unique(as.character(seqnames(zz)))]]
    score <- Views( cov, start = start(zz), end = end(zz) ) %>% sum()
    tibble(wig = my.wigs,value = score,rowname = zz$name)
  }) %>% bind_rows()
}
gamma.bw <- PhDfunc::GetBWList()[["GAM_clean_normalized_01022018_pOHT"]] %>% import.bw(as="RleList")
obserOE <- "oe"
# expe <- c("HiC_D_DIvA","HiC_D_OHT","HiC_D_OHTATMi","HiC_D_OHTDNAPKi","HiC_D_OHTPARPi")
expe <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_BENJ/data/",pattern=".hic") %>% str_remove(".hic")
# dir_dump <- "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/"
dir_dump <- "/home/rochevin/Documents/PROJET_INGE/HiC_BENJ/dump/TRANS/KR"
# Resolution
binSize=100000
binSizeTxt=paste0(binSize/1e3,"kb")

SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)

my_chromosomes <- seqlevels(DSB174)


mes_combi <- tibble(
  combi1 = expe[c(4,2)],
  combi2 = expe[c(3,1)]
) %>% as.matrix()
for(j in 1:nrow(mes_combi)){
  
  DIvA <- mes_combi[j,2]
  # OHT <- expe[4]
  OHT <- mes_combi[j,1]
  
  res_pca <- mclapply(my_chromosomes,function(Chr.V){
    message(Chr.V)
    chrendi=seqlengths(SeqInfo[Chr.V])
    binendi=ceiling(chrendi/binSize)
    
    fileCounti_DIvA=glue::glue("{dir_dump}/dump_{obserOE}_KR_{DIvA}_{Chr.V}_{Chr.V}_{binSizeTxt}.txt.gz")
    # fileCounti_OHT=paste0("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/data/KR/dump_",obserOE,"_KR_",expe[3],"_",Chr.V,"_",Chr.V,"_",binSizeTxt,".txt.gz")
    fileCounti_OHT=glue::glue("{dir_dump}/dump_{obserOE}_KR_{OHT}_{Chr.V}_{Chr.V}_{binSizeTxt}.txt.gz")
    
    
    
    data_mOHT <- read_Dumped_matrix(fileCounti_DIvA,binendi=binendi) %>% 
      as_tibble %>% 
      setNames(c("bin1","bin2","mOHT"))
    data_pOHT <- read_Dumped_matrix(fileCounti_OHT,binendi=binendi) %>% 
      as_tibble() %>%
      setNames(c("bin1","bin2","pOHT")) %>% inner_join(data_mOHT,by = c("bin1","bin2"))
    
    datai <- data_pOHT %>% mutate(ratio = log2(pOHT/mOHT))
    datai <- datai %>% mutate(ratio = ifelse(pOHT==0 & mOHT== 0,0,ratio))
    datai <- datai[is.finite(datai$ratio),]
    data.Mati=sparseMatrix(i=(datai$bin1/binSize)+1,j=(datai$bin2/binSize)+1,x=datai$ratio)
    data.MatSymi=data.Mati+t(data.Mati)
    diag(data.MatSymi)=diag(data.MatSymi)/2
    rm(data.Mati,datai)
    
    starti=seq(1,chrendi,by=binSize)
    endi=c(seq(binSize,chrendi,by=binSize),chrendi)
    chri.GR=GRanges(Chr.V,IRanges(starti,endi))
    names(chri.GR)=paste0("bin",1:length(chri.GR))
    
    HTC=HTCexp(data.MatSymi,chri.GR,chri.GR)
    
    cm <- HiTC:::sparseCor(intdata(HTC))
    cm[is.na(cm)] <- 0
    princp = princomp(cm)
    
    
    PCA_1 <-bind_cols(as_tibble(x_intervals(HTC)),as_tibble(loadings(princp)[])) %>%
      dplyr::select(seqnames:strand,str_c("Comp",1:10,sep = ".")) %>% 
      as_granges()%>% join_overlap_left(DSB174)  %>% sortSeqlevels() %>% sort()
    
    gamma.val <- Get1val(x = PCA_1,my.wigs = "gamma",one.w = gamma.bw) 
    
    PCA_cor_DSB <- PCA_1 %>% mutate(DSB = ifelse(!is.na(name),1,0)) %>% as_tibble() %>% dplyr::select(-seqnames:-strand,-name,-score) %>% 
      mutate(gamma = gamma.val$value)
    
    CompNames <- PCA_cor_DSB %>% names() %>% str_subset("Comp")
    PCA_cor_DSB <- lapply(CompNames,function(one_n){
      one_col <- PCA_cor_DSB[[one_n]]
      cc <- cor.test(one_col,PCA_cor_DSB$gamma)
      tibble(my_cor= cc$estimate,p.val = cc$p.value,name = one_n)
    }) %>% bind_rows()
    
    
    
    PCA_1 <- PCA_1 %>% as_tibble() %>%
      gather(key = PCs,value = value,-seqnames:-strand,-name:-score) %>% 
      left_join(PCA_cor_DSB,by = c("PCs"="name"))
    
    PCA_1 %>% filter(PCs %in% str_c("Comp",1:10,sep = ".")) 
    
  },mc.cores=4) %>% bind_rows()
  
  
  
  
  res_pca <- res_pca  %>% as_granges()
  
  # cc1 <- res_pca %>% filter_by_overlaps(DSB174 %>% anchor_center() %>% mutate(width = 10000)) %>% as_tibble() %>% group_by(seqnames) %>%
  #   summarise(meanval = mean(value)) %>% filter(meanval < 0) %>% pull(seqnames) %>% as.character()
  # 
  # 
  # res_pca <- res_pca %>% as_tibble() %>%split(.,.$seqnames)
  # res_pca <- res_pca[-1]
  # res_pca <- res_pca %>% map(function(x){
  #   if(unique(x$seqnames) %in% cc1){
  #     x %>% mutate(value = value * -1)
  #   }else{
  #     x
  #   }
  # })  
  
  # res_pca <- res_pca  %>% as_granges()
  
  order_chr <- res_pca %>% as_tibble() %>%filter(PCs %in% str_c("Comp",1,sep = ".")) %>%  group_by(seqnames) %>% summarise(corrmoi = mean(my_cor)) %>% arrange(desc(corrmoi)) %>% pull(seqnames)
  
  res_pca <- res_pca %>% as_tibble() %>%split(.,.$seqnames)
  # res_pca <- res_pca[-1]
  
  res_p <- lapply(res_pca[order_chr],function(one_data){
    one_data <- one_data %>% filter(PCs %in% str_c("Comp",1,sep = "."))%>% mutate(my_cor = round(my_cor,3))%>% mutate(p.val = format.pval(p.val,3))%>% unite("PCs",PCs,my_cor,sep = " cor: ") %>% unite("PCs",PCs,p.val,sep = " p.val: ") 
    one_data %>%  ggplot(aes(x=start,y=value)) + geom_line() +
      facet_wrap(~PCs,scales="free",ncol=1) +
      geom_vline(data = filter(one_data,!is.na(name)),aes(xintercept=start),linetype ="dashed") +ggtitle(unique(one_data$seqnames))
  })
  
  pdf(glue::glue("results/compD/{OHT}_{DIvA}_test_pca_cor_alamano_100kb_1_with_reverse_signal_{binSizeTxt}.pdf"),height=4,width=12)
  print(res_p)
  dev.off()
  
  ##Extract to bigwig 
  cc <- res_pca %>% map(filter,PCs %in% str_c("Comp",1,sep = ".")) %>% map(mutate,value=ifelse(my_cor < 0,value*-1,value)) %>% 
    bind_rows() %>% dplyr::select(seqnames:strand,value) %>% as_granges()
  # cc <- res_pca %>% bind_rows() %>% filter(PCs == "Comp.1") %>% dplyr::select(seqnames:strand,value) %>% as_granges()
  
  cc <- regioneR::filterChromosomes(cc,keep.chr=my_chromosomes)
  seqlengths(cc) <- seqlengths(SeqInfo)[names(seqlengths(cc))]
  cc.cov <- coverage(cc,weight="value")
  export.bw(cc.cov,glue::glue("results/compD/PC1_all_chr_log2ratio_{binSizeTxt}_{OHT}_{DIvA}.bw"))
  
}

