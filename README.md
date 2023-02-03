# CircadianClock
-------------------------------

This repository contains scripts that were used generate the figures for our paper (add final paper title)

### Data availability 
Raw and processed data are available at

------------------------------


### Script Overview

| Scripts           | Description                                                                                       | Figures               |
|-------------------|---------------------------------------------------------------------------------------------------|-----------------------|
|Bash_scripts/APA.sh                    | Script to generate data for APA plots                     |                          |
|Bash_scripts/ChIP_seq_script_SUN1SUN2.sh | Pipeline to process raw ChIP-seq data |                                                             |
|Bash_scripts/deeptools.sh  |  Code to generate count matrices from deeptools to be used with ggplot   |                                        |
|scriptsR/APA_plot.r | Function and code to plot APA from juicer tools |                                            Figs. 5F, S5G      |
|scriptsR/AverageProfile_Boxplot_Heatmap_DSB.r | Function and code for DSB profiles, boxplots, and heatmaps from deeptools counts matrices |  Figs 3C, 3D , 4E, S3B, S3E, S4H      |
|scriptsR/AverageProfile_GeneBodies.R| Function and code for average profiles at gene bodies |    Fig S3C, S4G    |

--------------------------------

### Data Overview

| Data        | Description                                                                                       |
|-------------------|---------------------------------------------------------------------------------------------------|
|data/80random.bed                 | 80 random sites                  |
|data/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed | 80 best AsiS1 DSB cut sites | 
|data/BLESS_HR_JunFragPE_Rmdups_pm500bp.bed  |  HR DSBs   |
|data/BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed |  NHEJ DSBs  |
|data/BMAL1DIVA_q0.1_peaks.narrowPeak  |  MACS2 Peak calling output from BMAL1 ChIP-seq  |
