## plot avgprof, heatmap, and boxplot from deeptools matrix


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


##############################################
#### Plotting


## Formatting info:

## The inputs to be used in the function are as follows:
## FILES: The matrices that are outputed from deeptools - the matrices need to have only one bw and one bed per matrix.
## vars: the sample names to be used, in the same order as the matrix files,
## regions: the peak region names to be used, in the same order as the matrix files

## all inputs need to be in a vector ( this is represented by the c() ) with each in inverted commas and separated by a comma, for example:
## FILES = c("file1", "file2")

## run function to read in data
data <- plot_group(
  FILES = c("/media/scollins/Sarah_scRNA/PER2_BMAL1/bwC_PER2_80_DSB_10kb.ma.gz",
            "/media/scollins/Sarah_scRNA/PER2_BMAL1/bwC_PER2_80_random_10kb.ma.gz"), 
  vars = c(rep(c("PER2"),1)),
  regions=c("80_DSB","80_random"),
  strand = "same"
)


## save data to individual variables for plotting
dat.heatmap <- data[[1]]
dat.plot <- data[[2]]
dat.boxplot <- data[[3]]

## Plotting variables...

plot_title = "ChIP-seq bwC PER2 at 80 DSB \n log2(+DSB/-DSB) " # title at top of plot
y_axis_title=  "ChIPseq  \n log2(+DSB/-DSB)" # y axis label for average profile


var_order = c("PER2") # order to show samples (needs to be exact sample as "vars" input in function)
regions_order = c("80_DSB","80_random") # Order to show regions if more than one used. If not, just name the one
hm_order_var = "PER2" # sample used to order heatmap rows (needs to be exact sample as one of the "vars" input in function)
basepairs = 10000 #Total number of basepairs that is used 
brewer_palette="RdBu" # brewer color palette used for heatmap
filename="/media/scollins/Sarah_scRNA/PER2_BMAL1/plots_for_gaelle/bwC_PER2_SUN1_10kb" # output file name title - do not add suffix 

size <- paste0((basepairs / 1000) / 2,"kb")

## MAKE HEATMAP

## order rows on all heatmaps based on selected variable
ord <- subset(dat.heatmap, variable==hm_order_var) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(value=mean(value))
ord <- ord$site[order(-ord$value)]

dat.heatmap$site <- factor(dat.heatmap$site, levels=rev(ord))


dat.heatmap$variable <- factor(dat.heatmap$variable, levels=c(var_order)) # set variable factors
dat.heatmap$Peaks <- factor(dat.heatmap$Peaks, levels=c(regions_order)) # set regions factors


hm_cutoff = quantile(dat.heatmap$value, 0.99)# heatmap cutoff value (sometimes you will need to scale the heatmaps properly by removing large outlier values)
# heatmap cutoff value (sometimes you will need to scale the heatmaps properly by removing large outlier values)


hm_cutoff_m  = quantile(dat.heatmap$value, 0.01)
# ## set heatmap cutoff if selected
if(!is.null(hm_cutoff)){
  dat.heatmap$value[dat.heatmap$value>hm_cutoff] <- hm_cutoff
  dat.heatmap$value[dat.heatmap$value<hm_cutoff_m] <- hm_cutoff_m
} else{
  dat.heatmap <- dat.heatmap
}

lim <- pmax(abs(hm_cutoff), abs(hm_cutoff_m))


p<- ggplot(dat.heatmap, aes(window, site)) +
  geom_raster(aes(fill = value, color=value)) + 
  facet_grid(Peaks~variable, scales = "free_y", space = "free_y") +
  ylab(" ") +
  scale_x_discrete(name = paste0("Distance from DSB (bp)"),
                   breaks = c(1, length(levels(as.factor(dat.heatmap$window))) / 2 ,length(levels(as.factor(dat.heatmap$window)))),
                   labels = c(paste0("-",size[1]), 'DSB', paste0("+",size[1]))) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_text(size=5),
        axis.text.x = element_text(size=8, face="bold"),
        axis.title.x = element_text(size=8),
        plot.title = element_text(size=7, face="bold", hjust=0.5),
        legend.position = "right") +
  scale_fill_distiller(palette = brewer_palette, direction = -1, na.value = "red", limits=c(-lim,lim)) +
  scale_color_distiller(palette = brewer_palette, direction = -1, na.value = "red", limits=c(-lim,lim)) + 
  ggtitle(plot_title) + facet_grid(Peaks~variable, scales="free_y")

print(p)

width_hm = 4.5 # width in inches for heatmap file
height_hm = 5  # height in inches for heatmap file

pdf(paste0(filename,"_heatmap.pdf"), width= width_hm, height=height_hm)
print(p)
dev.off()



### MAKE AVERAGE PROFILE

dat.plot$Peaks <- factor(dat.plot$Peaks,  levels=c(regions_order))
dat.plot$variable <- factor(dat.plot$variable, levels=c(var_order))

p <-ggplot( na.omit( dat.plot), aes( Window, Value) ) +
  labs( list( title = "", x = "", y = " " )) +
  geom_line(aes(color=variable), size=0.4, alpha=0.7) +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size=9),
        axis.text.x = element_text(size=9),
        axis.title.x = element_text(size=9),
        plot.title = element_text(size=9, face="bold", hjust=0.5),
        legend.position = "bottom") + ylab(y_axis_title) +
  
  scale_x_continuous(name = paste0("Distance from DSB (bp)"),
                     breaks = c(1, length(levels(as.factor(dat.plot$Window))) / 2 ,length(levels(as.factor(dat.plot$Window)))),
                     labels = c(paste0("-",size[1]), 'DSB', paste0("+",size[1]))) +
  geom_hline(yintercept = 0, linetype="dashed") + 
  geom_vline(xintercept = length(levels(as.factor(dat.plot$Window))) / 2, linetype="dashed", color="grey") +
  scale_color_manual(values=c("navy","firebrick4")) +
  ggtitle(paste0(plot_title))

print(p)

width_avgP = 3.5 # width in inches for average profile file
height_avgP = 3.5 # height in inches for average profile file

pdf(paste0(filename,"_AvgProf.pdf"),width=width_avgP, height=height_avgP)
print(p)
dev.off()


##### MAKE BOXPLOT

dat.boxplot$variable <- factor(dat.boxplot$variable, levels=c(var_order)) # set variable factors
dat.boxplot$Peaks <- factor(dat.boxplot$Peaks, levels=c(regions_order)) # set regions factors

stat.test <- as.data.frame(dat.boxplot) %>%
   wilcox_test(Value ~ Peaks*variable) %>%
   add_significance()
 stat.test

g <- ggplot(na.omit(dat.boxplot), aes(x = Peaks, y = (Value), fill=Peaks)) + 
  geom_boxplot() + 
  ggtitle("") + 
  theme_classic() +
  theme( 
    axis.text.x = element_text(size=11, angle=45, hjust=1),
    axis.text.y  = element_text(size=12),
    axis.title = element_text(angle = 0, hjust = 0.5, size=15),
    legend.position = "none",
    legend.text = element_text(size=13),
    legend.title = element_text(size=13),
    plot.title = element_text(hjust = 0.5, face="bold", size=15)) +
  ylab(paste0(y_axis_title)) + 
  facet_wrap(~variable) +
  xlab("Condition") 

print(g)

width_avgP = 6 # width in inches for average profile file
height_avgP = 5 # height in inches for average profile file

pdf(paste0(filename,"_boxplot.pdf"),width=width_avgP, height=height_avgP)
print(g)
dev.off()
