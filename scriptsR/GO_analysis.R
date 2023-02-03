require(tidyverse)
require(clusterProfiler)
require(org.Hs.eg.db)
require(ggtext)

theme_set(theme_minimal(base_size=12))
theme_update(
  legend.position = "bottom",axis.text.y = ggtext::element_markdown(),plot.title = ggtext::element_markdown()
)
files_csv <- list.files(pattern=".csv")

data_genes <- files_csv %>% setNames(str_remove(files_csv,".csv")) %>% map(read_tsv,col_names = c("Gene_ID")) %>% 
  bind_rows(.id = "File")


conversion_id <- bitr(data_genes$Gene_ID,fromType = "SYMBOL",toType = "ENTREZID",OrgDb=org.Hs.eg.db)

data_genes <- data_genes %>% left_join(conversion_id,by = c("Gene_ID"="SYMBOL")) %>% drop_na()



formula_res <- compareCluster(ENTREZID~File, data=data_genes, fun="enrichGO",ont="All",OrgDb=org.Hs.eg.db,
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05)


res_data <- fortify(formula_res, showCategory=NULL, by="Count",
                    includeAll=TRUE, split=NULL)
require(ggtext)

res_data <- res_data %>%
  arrange(Description) %>% mutate(Order = 1:dplyr::n()) %>% 
  mutate(color = case_when(
    ONTOLOGY == "BP" ~ "#27ae60",
    ONTOLOGY == "CC" ~ "#2980b9",
    ONTOLOGY == "MF" ~ "#d35400",
  )) %>% 
  mutate(
    name = glue::glue("<i style='color:{color}'>{Description}</i> ({ID})")
  )  %>% rowwise() %>% mutate(EvalGeneRatio = eval(parse(text =GeneRatio)))

theme_set(theme_minimal(base_size=12))
theme_update(
  legend.position = "bottom",axis.text.y = ggtext::element_markdown(),plot.title = ggtext::element_markdown()
)


res_data.splited <- split(res_data,res_data$File)

for(i in names(res_data.splited)){
  sres <- res_data.splited[[i]]
  p <- sres %>% ggplot(aes(x=Count,y=fct_reorder(name,Count,median),col=color,label = Count)) +
    geom_segment(aes(yend=fct_reorder(name,Order),xend=0,x=Count),linetype="dashed") +
    gggibbous::geom_moon(
      aes(ratio = EvalGeneRatio,fill=color),size=6,
      right = FALSE,
      key_glyph = gggibbous::draw_key_moon_left
    ) +
    gggibbous::geom_moon(
      aes(ratio = 1-EvalGeneRatio),fill="black",size=6
    ) +
    geom_text(size=3.5,col="white", fontface = "bold") +
    scale_color_identity() + 
    scale_fill_identity() 

  
  p <- p +
    labs(
      title = i,
      x = "Gene count", 
      y = "Ontologies "
    )
  outname <- str_replace_all(i,'\\s','_')
  ggsave(glue::glue("{outname}_GO_ontologies.pdf"),p,height=45,width=18)
}

##Take only BP and MF now
## All
data_genes.all <- data_genes  %>% filter(str_detect(File,"Banana",negate=T))



formula_res_all_BP <- compareCluster(ENTREZID~File, data=data_genes.all, fun="enrichGO",ont="BP",OrgDb=org.Hs.eg.db,
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01)

formula_res_all_CC <- compareCluster(ENTREZID~File, data=data_genes.all, fun="enrichGO",ont="CC",OrgDb=org.Hs.eg.db,
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.01)



for(i in c(10)){
  res_data <- lapply(list("BP" = formula_res_all_BP,"CC"=formula_res_all_CC),function(x){
    
    x %>% dropGO(level = i) %>% simplify(cutoff=0.7, by="p.adjust", select_fun=min) %>% fortify(showCategory=NULL, by="Count",
                                                        includeAll=TRUE, split=NULL) 
    # bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
  }) %>% bind_rows(.id = "ONTOLOGY")
  
  
  res_data <- res_data %>%
    arrange(Description) %>% mutate(Order = 1:dplyr::n()) %>% 
    mutate(color = case_when(
      ONTOLOGY == "BP" ~ "#27ae60",
      ONTOLOGY == "CC" ~ "#2980b9",
      ONTOLOGY == "MF" ~ "#d35400",
    )) %>% 
    mutate(
      name = glue::glue("<i style='color:{color}'>{Description}</i> ({ID})")
    )  %>% rowwise() %>% mutate(EvalGeneRatio = eval(parse(text =GeneRatio)))
  
  
  p1 <- filter(res_data,File == "All targeted genes for GO") %>% ggplot(aes(x="All targeted genes for GO",y=fct_reorder(name,Count,median),color=p.adjust,label = Count)) +

    geom_point(aes(size=Count)) +

    scale_color_gradientn(colors=MetBrewer::met.brewer("Tam"),
      guide=guide_colorbar(title="Adjusted P-value",
                           title.position = "top",
                           title.hjust = .5,
                           barwidth = unit(10, "lines"),
                           barheight = unit(1, "lines")),breaks = c(0,0.005,0.01),limits = c(0,0.01)) +
    facet_wrap(~ONTOLOGY,scales="free",ncol=2)+
    labs(
      title = "All targeted genes for GO",
      x = "Gene count", 
      y = "Ontologies "
    ) + guides(size = guide_legend(title = "Gene ratio"))
  
  p2 <- filter(res_data,File != "All targeted genes for GO") %>%
    mutate(File = str_extract(File,"Decreasing|Increasing")) %>% 
    ggplot(aes(x=File,y=fct_reorder(name,Count,median),color=p.adjust,label = Count)) +
    
    geom_point(aes(size=Count)) +
    
    scale_color_gradientn(colors=MetBrewer::met.brewer("Tam"),
                          guide=guide_colorbar(title="Adjusted P-value",
                                               title.position = "top",
                                               title.hjust = .5,
                                               barwidth = unit(10, "lines"),
                                               barheight = unit(1, "lines")),breaks = c(0,0.005,0.01),limits = c(0,0.01)) +
    facet_wrap(~ONTOLOGY,scales="free",ncol=2)+
    labs(
      title = "All targeted genes for GO",
      x = "Gene count", 
      y = "Ontologies "
    ) + guides(size = guide_legend(title = "Gene ratio"))
  ggsave(glue::glue("AllTargeted_Golevel_{i}_ontologies.pdf"),p1,height=12,width=12)
  ggsave(glue::glue("GOMeanArea_Golevel_{i}_ontologies.pdf"),p2,height=10,width=16)
}



p <- p 
outname <- str_replace_all(i,'\\s','_')
ggsave(glue::glue("{outname}_GO_ontologies.pdf"),p,height=45,width=18)