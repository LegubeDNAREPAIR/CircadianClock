### Average Profile on Gene Bodies


require("dplyr")
require("ggplot2")
require("reshape2")
require("rtracklayer")
library("BSgenome.Hsapiens.UCSC.hg19")
library("plyranges")
library("tidyverse")
library("RColorBrewer")
seqlens = seqlengths( Hsapiens );

##################
## Functions for taking deeptools matrix and making compatible for ggplot
#################


plot_group <- function(FILES, vars, regions, strand="same"){
  
  dat.plot <- NULL
  dat.heatmap <- NULL
  dat.boxplot <- NULL
  
  
  message("Reading in data...")
  ## load in data to dfs
  for(n in 1:length(FILES)){
    df <- read.csv(FILES[n], sep="\t", header = FALSE, skip=1)
    rownames(df) <- df$V4
    #rownames(df) <- paste0(df$V1,":",df$V2,"-",df$V3)
    df <- (df[,-c(1:6)])
    df[df=="NaN"] <- 0
    
    if(strand=="split"){
      df <- df[,ncol(df):1]
    }else{
      df <- df
    }
    
    if(nrow(df)>400){
      s <- nrow(df) / 400
      df <- df[ order(rowMeans(df)), ]
      df <- aggregate(df, list(rep(1:(nrow(df) %/% s + 1), each = s, len = nrow(df))), mean)[-1];
    }
    
    sub.dat.plot <- data.frame(
      Window=c(1:ncol(df)),
      Value=colMeans((df)),
      variable = vars[n],
      Peaks=regions[n]
    )
    
    if(is.null(dat.plot)){
      dat.plot <- sub.dat.plot
    }else{
      dat.plot <- rbind(dat.plot,sub.dat.plot)
    }
    
    
    colnames(df) <- 1:ncol(df)
    melted <- melt(df)
    colnames(melted) <- c("window","value")
    melted$variable <- vars[n]
    melted$Peaks <- regions[n]
    melted$site <- rownames(df)
    
    if(is.null(dat.heatmap)){
      dat.heatmap <- melted
    }else{
      dat.heatmap <- rbind(dat.heatmap,melted)
    }   
    
    
    sub.dat.boxplot <- data.frame(
      Value=rowMeans(df),
      variable = vars[n],
      Peaks=regions[n]
    )
    
    if(is.null(dat.plot)){
      dat.boxplot <- sub.dat.boxplot
    }else{
      dat.boxplot <- rbind(dat.boxplot,sub.dat.boxplot)
    }
    
  }
  
  return(list(dat.heatmap,dat.plot,dat.boxplot))
  
} 

###########################################################################################

for(cond in c("SUN2","SUN1")){
  for(bp in c(2000)){
    
    sizekb <- paste0((bp / 1000),"kb")
    regions="all_ensembl_genes"
    
    ## run function to read in data
    data <- plot_group(
      FILES = c(paste0("/media/scollins/Sarah_scRNA/PER2_BMAL1/",cond,"DIVA_EnsDB_genes_plus_",sizekb,".ma.gz"),
                paste0("/media/scollins/Sarah_scRNA/PER2_BMAL1/",cond,"OHT_EnsDB_genes_plus_",sizekb,".ma.gz")
                
      ), 
      vars = c(paste0(cond,"_DIVA"),
               paste0(cond,"_OHT")
      ),
      regions=c(rep(c("allEnsDB_genes"),2)
      ),
      strand = "same"
    )
    
    
    # c(rep("80_DSB",2))
    
    
    ## save data to individual variables for plotting
    dat.heatmap <- data[[1]]
    dat.plot <- data[[2]]
    dat.boxplot <- data[[3]]
    
    data <- plot_group(
      FILES = c(paste0("/media/scollins/Sarah_scRNA/PER2_BMAL1/",cond,"DIVA_EnsDB_genes_minus_",sizekb,".ma.gz"),
                paste0("/media/scollins/Sarah_scRNA/PER2_BMAL1/",cond,"OHT_EnsDB_genes_minus_",sizekb,".ma.gz")
                
      ), 
      vars = c(paste0(cond,"_DIVA"),
               paste0(cond,"_OHT")
      ),
      regions=c(rep(c("allEnsDB_genes"),2)
      ),
      strand = "split"
    )
    
    dat.heatmap_r <- data[[1]]
    dat.plot_r <- data[[2]]
    dat.boxplot_r <- data[[3]]
    
    dat.heatmap <- rbind(dat.heatmap,dat.heatmap_r)
    dat.plot <- rbind(dat.plot,dat.plot_r)
    dat.boxplot <- rbind(dat.boxplot,dat.boxplot_r)
    
    
    dat.plot <- dat.plot %>% 
      dplyr::group_by(Window,variable,Peaks) %>% 
      dplyr::summarise(Value=mean(Value))
    
    
    dat.plot <- (dat.plot) %>% tidyr::separate(variable,c("variable","trt"), sep="_(?=[^_]+$)")
    dat.heatmap <- (dat.heatmap) %>% tidyr::separate(variable,c("variable","trt"), sep="_(?=[^_]+$)")
    dat.boxplot <- (dat.boxplot) %>% tidyr::separate(variable,c("variable","trt"), sep="_(?=[^_]+$)")
    
    ## Plotting variables...
    
    plot_title = paste0("ChIP-seq ",cond," ", sizekb," at ",regions) # title at top of plot
    y_axis_title=  paste0("ChIP-seq  ",cond," \n Avg Norm Read Count") # y axis label for average profile
    
    var_order = c(cond) # order to show samples (needs to be exact sample as "vars" input in function)
    regions_order = c("allEnsDB_genes") # Order to show regions if more than one used. If not, just name the one
    hm_order_var = cond # sample used to order heatmap rows (needs to be exact sample as one of the "vars" input in function)
    basepairs = bp #Total number of basepairs that is used 
    brewer_palette="RdBu" # brewer color palette used for heatmap
    filename=paste0("/media/scollins/Sarah_scRNA/PER2_BMAL1/",cond,"_",regions,"_",sizekb) # output file name title - do not add suffix 
    
    size <- paste0((bp / 1000) / 2,"kb")
    
    dat.plot$Peaks <- factor(dat.plot$Peaks,  levels=c(regions_order))
    dat.plot$variable <- factor(dat.plot$variable, levels=c(var_order))
    dat.plot$trt <- factor(dat.plot$trt, levels=c("DIVA","OHT"))
    
    mycolor <- c("#2980b9","#c0392b")
    
    avgp <-ggplot( na.omit( dat.plot), aes( Window, Value, color=trt) ) +
      labs( list( title = "", x = "", y = " " )) +
      geom_line(size=0.6, alpha=0.7) +
      theme_classic(base_size = 10) +
      theme(axis.text.y = element_text(size=9),
            axis.text.x = element_text(size=9),
            axis.title.x = element_text(size=7),
            plot.title = element_text(size=9, face="bold", hjust=0.5),
            legend.position = "bottom") + ylab(y_axis_title) +
      
      scale_x_continuous(name = paste0("Genes TSS -> TES (bp)"),
                         breaks = c(1, 40 , 140, 180),
                         labels = c("-2kb","TSS","TES","+2kb")
      ) +
      geom_vline(xintercept = 40, linetype="dashed", color="grey") +
      geom_vline(xintercept = 140, linetype="dashed", color="grey") +
      scale_color_manual(values = mycolor,name="Pol2 Cat") +
      ggtitle(paste0(plot_title))
    
    print(avgp)
    
    pdf(paste0(filename,"_AvgProf.pdf"),width=6, height=3.5)
    print(avgp)
    dev.off()
    
  
    
  }
}

