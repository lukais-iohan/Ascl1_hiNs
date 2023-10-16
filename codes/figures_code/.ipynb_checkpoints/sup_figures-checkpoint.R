#### Supp figures

library("Seurat")
library("SingleCellExperiment")
library("plyr")
library('dplyr')
library("DESeq2")
library('ggplot2')
library("gprofiler2")
library('scales')
library('ggrepel')
library('ggplot2')

EXC_DEGs <- read.csv("../data/EXC_DEGs.csv",header = T, stringsAsFactors = F, row.names = 1)


INH_DEGs <- read.csv("../data/INH_DEGs.csv",header = T, stringsAsFactors = F, row.names = 1)


EXC_DEGs <- EXC_DEGs %>% mutate(diff_expressed = case_when(
  avg_log2FC > 0.25 & p_val_adj < 0.01 ~ 'UP',
  avg_log2FC < -0.25 & p_val_adj < 0.01 ~ 'DOWN',
  TRUE ~ 'NO'
))


INH_DEGs <- INH_DEGs %>% mutate(diff_expressed = case_when(
  avg_log2FC > 0.25 & p_val_adj < 0.01 ~ 'UP',
  avg_log2FC < -0.25 & p_val_adj < 0.01 ~ 'DOWN',
  TRUE ~ 'NO'
))


EXC_DEGs %>% filter(comparasion %in% c('hiNPCs_6w x hiNPCs_4w','hiNPCs_6w x NEUROG2','NEUROG2 x hiNPCs_4w')) %>%
  ggplot(aes(x=avg_log2FC, y=-log10(p_val_adj), col=diff_expressed)) +
  geom_point(aes(size = 1.5),alpha = 0.5,show.legend = F) + 
  theme_ipsum(base_size = 20,axis_text_size = 35,plot_title_size = 50,grid = F,axis_title_size = 20,axis_title_family = 20)+
  scale_color_manual(values=c("dodgerblue2", "black", "firebrick2")) +
  geom_vline(xintercept=c(-0.25, 0.25), col="black",linetype = 'dashed') +
  ylim(0,100)+ xlim(-2.5,2.5) + xlab('') + ylab('')+
  geom_hline(yintercept=-log10(0.01), col="black",linetype = 'dashed') +facet_wrap(~comparasion,nrow = 3,ncol = 1,scales = 'free') 


ggsave('figuras_artigo/fig2/EXC/exc_degs_sup.jpeg',width = 12,height = 12,dpi = 600)


INH_DEGs %>% filter(comparasion %in% c('hiNPCs_6w x hiNPCs_4w','hiNPCs_6w x NEUROG2','NEUROG2 x hiNPCs_4w')) %>%
  ggplot(aes(x=avg_log2FC, y=-log10(p_val_adj), col=diff_expressed)) +
  geom_point(aes(size = 1.5),alpha = 0.5,show.legend = F) + 
  theme_ipsum(base_size = 20,axis_text_size = 35,plot_title_size = 50,grid = F,axis_title_size = 20,axis_title_family = 20)+
  scale_color_manual(values=c("dodgerblue2", "black", "firebrick2")) +
  geom_vline(xintercept=c(-0.25, 0.25), col="black",linetype = 'dashed') +
  ylim(0,100)+ xlim(-2.5,2.5) +xlab('') + ylab('')+
  geom_hline(yintercept=-log10(0.01), col="black",linetype = 'dashed') +facet_wrap(~comparasion,nrow = 3,ncol = 1,scales = 'free') 

ggsave('figuras_artigo/fig2/INH/inh_degs_sup.jpeg',width = 10,height = 12,dpi = 600)



### GSEA Analysis
GSEA = function(gene_list, GO_file, pval,comparasion) {
  library(dplyr)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  
  fgRes <- fgsea::fgseaMultilevel(pathways = GO_file,
                                  stats = gene_list,
                                  minSize=15, ## minimum gene set size
                                  maxSize=400, ## maximum gene set size
                                  nPermSimple = 10000
  ) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = GO_file,
                                      stats = gene_list)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  
  
  if( nrow(fgRes) >0 ) {
    
    fgRes$'comparasion' <- comparasion
    
    return(fgRes)
  }
  
}

library('Seurat')
library('msigdbr')
library('singleseqgset')
library('heatmap3')
library('CEMiTool')
library('tidyr')
library('dplyr')
library('stringr')
library('ggplot2')

gmt_in <- read_gmt('data/c5.all.v2023.1.Hs.symbols.gmt') ### arquivo enviado em anexo
gmt_in <- gmt_in[str_detect(gmt_in$term,pattern = 'GOBP|GOCC|GOMF'),]
gmt_in$term <- str_remove_all(gmt_in$term,'GOBP|GOCC|GOMF')
gmt_in$term <- str_replace_all(gmt_in$term,'_',' ')
gmt_in <- split(gmt_in$gene,gmt_in$term)

dataset_degs_list <- list(EXC_DEGs,INH_DEGs)

gsea_datasets <- lapply(dataset_degs_list, function(x) {
  
  gsea_function <- function(x,comparasion) { ### function to get genes and sort them by avg_log2FC
    
    data_frame <- x[x$comparasion == comparasion,] 
    gene <- data_frame$metric
    names(gene) <- data_frame$genes
    gene <- sort(gene,decreasing = T)
    
    
  }
  
  comparasion_vector <- unique(x$comparasion)
  
  gene_list <- lapply(seq_along(comparasion_vector),function(i) { #### here we are runnning the gsea_function on ever comparasion
    ### the result of this code is a list with a list of genes in each comparasion
    
    res = gsea_function(x,comparasion_vector[i])
    
  })
  
  
  names(gene_list) <- comparasion_vector #### naming the comparasion. Note that the order of comparasion is the same you used in the comparasion vector
  
  
  ### Finnaly we run the GSEA analysis
  
  
  
  results_fgsea <- lapply(seq_along(comparasion_vector), function(i) { ### runnning the GSEA function with the genes list of comparasions
    
    res = GSEA(gene_list[[i]], gmt_in, pval = 0.01,comparasion_vector[i])
    
  })
  
  results_fgsea_sig <- results_fgsea[!sapply(results_fgsea,is.null)] 
  
  results_fgsea_sig <- Reduce(full_join,results_fgsea_sig)
  
})  



plot_function <- function(dataset,comparasion_dataset,neuron) {
  
  results <- dataset[dataset$comparasion == comparasion_dataset,]
  results_top <- bind_rows(
    
    results %>% group_by(comparasion) %>%
      slice_max(NES,n = 10),
    
    results %>% group_by(comparasion) %>%
      slice_min(NES,n = 10)
  )
  
  colos = setNames(c("firebrick2", "dodgerblue2"),
                   c("Up-regulated", "Down-regulated"))
  a <- ggplot(results_top,aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip()+ xlim(-3,3)
    labs(x="Pathway", y="Normalized Enrichment Score") + theme_bw(base_size = 10) +
    theme_bw()+ ggtitle(paste(comparasion_dataset,neuron,sep = " ")) +
    theme(axis.text.y = element_text(size = 12,face = 'bold'),
          axis.text.x = element_text(size = 12,face = 'bold'),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 12)) 
  return(a)
  
}

fgsea_EXC <- gsea_datasets[[1]]
fgsea_INH <- gsea_datasets[[2]]

### EXC

comparasion_vector <- unique(EXC_DEGs$comparasion)

lapply(seq_along(comparasion_vector), function(i) {
  
  p1 <- plot_function(fgsea_EXC,comparasion_vector[i],'EXC')
  ggsave(paste0('figuras_artigo/fig2/EXC/plot_',i,'.jpeg'),plot = p1,width = 12,height = 5,dpi = 600)
  
})

### INH

lapply(seq_along(comparasion_vector), function(i) {
  
  p1 <- plot_function(fgsea_INH,comparasion_vector[i],'INH')
  ggsave(paste0('figuras_artigo/fig2/INH/plot_',i,'.jpeg'),plot = p1,width = 15,height = 5,dpi = 600)
  
})


#### Intersection

library('genekitr')

### EXC

EXC_DEGs_diff <- EXC_DEGs %>% filter(diff_expressed != 'NO')

EXC_DEGs_list <- list('hiNPCs_6w x hiNPCs_4w' = EXC_DEGs_diff[EXC_DEGs_diff$comparasion == 'hiNPCs_6w x hiNPCs_4w',7],
                      'hiNPCs_6w x NEUROG2' = EXC_DEGs_diff[EXC_DEGs_diff$comparasion == 'hiNPCs_6w x NEUROG2',7],
                      'NEUROG2 x hiNPCs_4w' = EXC_DEGs_diff[EXC_DEGs_diff$comparasion == 'NEUROG2 x hiNPCs_4w',7])

jpeg('figuras_artigo/fig2/EXC/intersection_sup_EXC.jpeg',width = 12,height = 8,res = 600,units = 'in')

upset(fromList(EXC_DEGs_list),order.by = 'freq',set_size.numbers_size = T,text.scale = c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5),point.size = 5,
      sets.bar.color = c('#FF8080FF','#FF2A2AFF','#00FF00FF'))

dev.off()

### INH

INH_DEGs_diff <- INH_DEGs %>% filter(diff_expressed != 'NO')

INH_DEGs_list <- list('hiNPCs_6w x hiNPCs_4w' = INH_DEGs_diff[INH_DEGs_diff$comparasion == 'hiNPCs_6w x hiNPCs_4w',7],
                      'hiNPCs_6w x NEUROG2' = INH_DEGs_diff[INH_DEGs_diff$comparasion == 'hiNPCs_6w x NEUROG2',7],
                      'NEUROG2 x hiNPCs_4w' = INH_DEGs_diff[INH_DEGs_diff$comparasion == 'NEUROG2 x hiNPCs_4w',7])

jpeg('figuras_artigo/fig2/INH/intersection_sup_INH.jpeg',width = 12,height = 8,res = 600,units = 'in')

upset(fromList(INH_DEGs_list),order.by = 'freq',set_size.numbers_size = T,text.scale = c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5),point.size = 5,
      sets.bar.color = c('#FF8080FF','#FF2A2AFF','#00FF00FF'))

dev.off()
