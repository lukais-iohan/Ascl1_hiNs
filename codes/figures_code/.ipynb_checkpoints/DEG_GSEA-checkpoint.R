library('ggrepel')
library('tidyr')
library('qvalue')
library('hrbrthemes')
library('dplyr')
library('RColorBrewer')
library('ggplot2')
library('stringr')

#### Import of DEGs table



EXC_DEGs <- read.csv("only_neurons_tables/EXC_DEGs.csv",header = T, stringsAsFactors = F, row.names = 1)


INH_DEGs <- read.csv("only_neurons_tables/INH_DEGs.csv",header = T, stringsAsFactors = F, row.names = 1)


### Defining differential expressed genes by log2 fold change and adjusted p_value

EXC_DEGs <- EXC_DEGs %>% mutate(diff_expressed = case_when(
  avg_log2FC > 0.25 & p_val_adj < 0.01 ~ 'UP',
  avg_log2FC < -0.25 & p_val_adj < 0.01 ~ 'DOWN',
  TRUE ~ 'NO'
))

EXC_DEGs$'Neuro_Type' <- 'Excitatory'


INH_DEGs <- INH_DEGs %>% mutate(diff_expressed = case_when(
  avg_log2FC > 0.25 & p_val_adj < 0.01 ~ 'UP',
  avg_log2FC < -0.25 & p_val_adj < 0.01 ~ 'DOWN',
  TRUE ~ 'NO'
))

INH_DEGs$'Neuro_Type' <- 'Inhibitory'

### save only degs 
write.csv(EXC_DEGs,"../data/EXC.csv")
write.csv(INH_DEGs,"../data/INH.csv")


### Plot volcano plot 
EXC_DEGs$'label' <- EXC_DEGs$genes %in% c('SYN3',
                                          'SCN7A',
                                          'ANK1',
                                          'CACNA1A',
                                          'KCNQ1',
                                          'KCNQ2',
                                          'DLG2',
                                          'GRIA2',
                                          'GABRB1',
                                          'CAMK4',
                                          'SCN1A',
                                          'SCN2A',
                                          'APP',
                                          'COX8A',
                                          'COX7C',
                                          'COX6A1',
                                          'GRIA1',
                                          'NDUFA4',
                                          'DSCAM') ### Channel, synapse genes 
EXC_DEGs <- EXC_DEGs %>% mutate(label = case_when(
  label == 'TRUE' ~ genes,
  TRUE ~ ''
))

EXC_DEGs %>% filter(comparasion %in% c('ASCL1 x hiNPCs_4w','ASCL1 x hiNPCs_6w','ASCL1 x NEUROG2')) %>%
ggplot(aes(x=avg_log2FC, y=-log10(p_val_adj),label=label, col=diff_expressed)) +
  geom_point(aes(size = 0.2),alpha = 0.5,show.legend = F) + 
  geom_text_repel(box.padding = 0.8, max.overlaps = Inf,color = 'black')+
  theme_ipsum(base_size = 20,axis_text_size = 35,plot_title_size = 50,grid = F,
              axis_title_size = 20,axis_title_family = 20)+
  scale_color_manual(values=c("dodgerblue2", "black", "firebrick2")) +
  geom_vline(xintercept=c(-0.25, 0.25), col="black",linetype = 'dashed') +
  ylim(0,100)+ xlim(-2.5,2.5) + xlab('') + ylab('')+
  geom_hline(yintercept=-log10(0.01), col="black",linetype = 'dashed') + 
  theme(axis.text = element_text(color = 'black'))+
 facet_wrap(~comparasion,nrow = 3,ncol = 1,scales = 'free') 


ggsave('figuras_artigo/fig2/EXC/exc_degs.jpeg',width = 10,height = 12,dpi = 600)

#### 


INH_DEGs$'label' <- INH_DEGs$genes %in% c('SYN3',
                                          'SCN7A',
                                          'ANK1',
                                          'CACNA1A',
                                          'KCNQ1',
                                          'KCNQ2',
                                          'DLG2',
                                          'GRIA2',
                                          'GABRB1',
                                          'CAMK4',
                                          'SCN1A',
                                          'SCN2A',
                                          'APP',
                                          'COX8A',
                                          'COX7C',
                                          'COX6A1',
                                          'GRIA1',
                                          'NDUFA4',
                                          'DSCAM')
INH_DEGs <- INH_DEGs %>% mutate(label = case_when(
  label == 'TRUE' ~ genes,
  TRUE ~ ''
))

INH_DEGs %>% filter(comparasion %in% c('ASCL1 x hiNPCs_4w','ASCL1 x hiNPCs_6w','ASCL1 x NEUROG2')) %>%
  ggplot(aes(x=avg_log2FC, y=-log10(p_val_adj),label=label, col=diff_expressed)) +
  geom_point(aes(size = 0.2),alpha = 0.5,show.legend = F) + 
  geom_text_repel(box.padding = 1, max.overlaps = Inf,color = 'black',force = 0.5)+
  theme_ipsum(base_size = 20,axis_text_size = 35,plot_title_size = 50,grid = F,
              axis_title_size = 20,axis_title_family = 20)+
  scale_color_manual(values=c("dodgerblue2", "black", "firebrick2")) +
  geom_vline(xintercept=c(-0.25, 0.25), col="black",linetype = 'dashed') +
  ylim(0,100)+ xlim(-2.5,2.5) + xlab('') + ylab('')+
  geom_hline(yintercept=-log10(0.01), col="black",linetype = 'dashed') + 
  theme(axis.text = element_text(color = 'black'))+
  facet_wrap(~comparasion,nrow = 3,ncol = 1,scales = 'free') 

ggsave('figuras_artigo/fig2/INH/inh_degs.jpeg',width = 10,height = 12,dpi = 600)




#### GSEA Analysis
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
library('CEMiTool')
library('tidyr')
library('stringr')
library('ggplot2')

gmt_in <- read_gmt('data/c5.all.v2023.1.Hs.symbols.gmt') ### gmt ontologie files 
gmt_in <- gmt_in[str_detect(gmt_in$term,pattern = 'GOBP|GOCC|GOMF'),] ### get only Gene Ontologies 
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

    
### Function to plot ontologies based on NES score 

plot_function <- function(dataset,comparasion_dataset) {
  
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
    coord_flip()+
    labs(x="Pathway", y="Normalized Enrichment Score") + theme_bw(base_size = 10) +
    theme_bw()+ ggtitle(comparasion_dataset) +
    theme(axis.text.y = element_text(size = 12,color = 'black'),
          axis.text.x = element_text(size = 12,color = 'black'),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 12),title = element_text(face = 'bold',colour = 'black')) 
  return(a)
  
}

fgsea_EXC <- gsea_datasets[[1]]

fgsea_EXC <- apply(fgsea_EXC, 2, as.character)
write.csv(fgsea_EXC,'only_neurons_tables/fgsea_EXC.csv')


fgsea_INH <- gsea_datasets[[2]]
fgsea_INH <- apply(fgsea_INH, 2, as.character)
write.csv(fgsea_INH,'only_neurons_tables/fgsea_INH.csv')

### EXC plot

comparasion_vector <- unique(EXC_DEGs$comparasion)

fgsea_plot_EXC <-lapply(seq_along(comparasion_vector), function(i) {
  
   plot_function(fgsea_EXC,comparasion_vector[i])
  
  
})


fgsea_plot_EXC[[3]]/fgsea_plot_EXC[[1]]/fgsea_plot_EXC[[2]]

ggsave('figuras_artigo/fig2/EXC/plot_3EXC.jpeg',width = 10,height = 15,dpi = 600)

fgsea_plot_EXC[[5]]/fgsea_plot_EXC[[4]]

ggsave('figuras_artigo/fig2/EXC/plot_sup_EXC.jpeg',width = 12,height = 13,dpi = 600)

## INH plot 

fgsea_plot_INH <-lapply(seq_along(comparasion_vector), function(i) {
  
  plot_function(fgsea_INH,comparasion_vector[i])
  
})

fgsea_plot_INH[[3]]/fgsea_plot_INH[[1]]/fgsea_plot_INH[[2]]

ggsave('figuras_artigo/fig2/INH/plot_3INH.jpeg',width = 10,height = 15,dpi = 600)

fgsea_plot_INH[[5]]/fgsea_plot_INH[[4]]/fgsea_plot_INH[[6]]

ggsave('figuras_artigo/fig2/INH/plot_sup_INH.jpeg',width = 10,height = 15,dpi = 600)

#### Intersection of genes 

library('UpSetR')

### EXC

EXC_DEGs_diff <- EXC_DEGs %>% filter(diff_expressed != 'NO')

EXC_DEGs_list <- list('ASCL1 x NEUROG2' = EXC_DEGs_diff[EXC_DEGs_diff$comparasion == 'ASCL1 x NEUROG2' ,7],
                      'ASCL1 x hiNPCs_4w' = EXC_DEGs_diff[EXC_DEGs_diff$comparasion == 'ASCL1 x hiNPCs_4w',7],
                      'ASCL1 x hiNPCs_6w' = EXC_DEGs_diff[EXC_DEGs_diff$comparasion == 'ASCL1 x hiNPCs_6w',7])


jpeg('figuras_artigo/fig2/EXC/intersection_EXC.jpeg',width = 12,height = 8,res = 600,units = 'in')

upset(fromList(EXC_DEGs_list),order.by = 'freq',set_size.numbers_size = T,text.scale = c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5),
      sets = c('ASCL1 x hiNPCs_4w','ASCL1 x hiNPCs_6w','ASCL1 x NEUROG2'),point.size = 5,sets.bar.color = c('#FF8080FF','#FF2A2AFF','#00FF00FF'))



dev.off()

### INH

INH_DEGs_diff <- INH_DEGs %>% filter(diff_expressed != 'NO')

INH_DEGs_list <- list('ASCL1 x NEUROG2' = INH_DEGs_diff[INH_DEGs_diff$comparasion == 'ASCL1 x NEUROG2' ,7],
                      'ASCL1 x hiNPCs_4w' = INH_DEGs_diff[INH_DEGs_diff$comparasion == 'ASCL1 x hiNPCs_4w',7],
                      'ASCL1 x hiNPCs_6w' = INH_DEGs_diff[INH_DEGs_diff$comparasion == 'ASCL1 x hiNPCs_6w',7])

jpeg('figuras_artigo/fig2/INH/intersection_INH.jpeg',width = 12,height = 8,res = 600,units = 'in')

upset(fromList(INH_DEGs_list),order.by = 'freq',set_size.numbers_size = T,
      text.scale = c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5),
      sets = c('ASCL1 x hiNPCs_4w','ASCL1 x hiNPCs_6w','ASCL1 x NEUROG2'),
      point.size = 5,sets.bar.color = c('#FF8080FF','#FF2A2AFF','#00FF00FF'))

dev.off()


##### Alternative way to plot ontologies

library(RedeR)
library(igraph)
library(TreeAndLeaf)

fgsea_EXC <- gsea_datasets[[1]]
fgsea_EXC_ASCL1 <- fgsea_EXC[fgsea_EXC$comparasion == 'ASCL1 x NEUROG2',] ### choose comparasion

gs_results <- fgsea_EXC_ASCL1
myGO <- gmt_in
df <- matrix(nrow=nrow(gs_results), ncol = nrow(gs_results), data = 0)
rownames(df) = colnames(df) = gs_results$pathway
  
for ( i in 1:nrow(gs_results)) {
    genesI =  unlist(myGO[names(myGO) == gs_results$pathway[i] ])
  for (j in 1:nrow(gs_results)) {
      genesJ = unlist(myGO[names(myGO) == gs_results$pathway[j] ])
      ## Jaccards distance  1 - (intersection / union )
      overlap = sum(!is.na(match(genesI, genesJ )))
      jaccards = overlap / length(unique(c(genesI, genesJ) ))
      df[i,j] = 1-jaccards
    }
  }
  
  ## Cluster nodes using dynamic tree cut, for colors
distMat = as.dist(df)
dendro = hclust(distMat, method = "average" )
clust = dynamicTreeCut::cutreeDynamicTree( dendro, minModuleSize = 4 )
  ## Note: in dynamicTreeCut, cluster 0, is a garbage cluster for things that dont cluster, so we remove it
  
gs_results$Cluster = clust

distMat = as.dist(df)
dendro = hclust(distMat, method = "average" )
  
tal <- treeAndLeaf(dendro)

tal <- att.mapv(g = tal, dat = gs_results, refcol = 1)

tal <- att.setv(g = tal, from = "Cluster", to = "nodeColor", 
                cols = c('blue','red'))
tal <- att.setv(g = tal, from = "size", to = "nodeSize")

#--- Set graph attributes using 'att.addv' and 'att.adde' functions
tal <- att.addv(tal, "nodeFontSize", value = 15, index = V(tal)$isLeaf)
tal <- att.adde(tal, "edgeWidth", value = 3)

rdp <- RedPort()
calld(rdp)
resetd(rdp)
#--- Send the tree-and-leaf to the interactive R/Java interface
addGraph(obj = rdp, g = tal, gzoom=50)



addLegend.color(obj = rdp, tal, title = "Group", 
                position = "topright",vertical=T)


####

fgsea_INH <- gsea_datasets[[2]]
fgsea_INH_ASCL1 <- fgsea_INH[fgsea_INH$comparasion == 'ASCL1 x NEUROG2',]

gs_results <- fgsea_INH_ASCL1
myGO <- gmt_in
df <- matrix(nrow=nrow(gs_results), ncol = nrow(gs_results), data = 0)
rownames(df) = colnames(df) = gs_results$pathway

for ( i in 1:nrow(gs_results)) {
  genesI =  unlist(myGO[names(myGO) == gs_results$pathway[i] ])
  for (j in 1:nrow(gs_results)) {
    genesJ = unlist(myGO[names(myGO) == gs_results$pathway[j] ])
    ## Jaccards distance  1 - (intersection / union )
    overlap = sum(!is.na(match(genesI, genesJ )))
    jaccards = overlap / length(unique(c(genesI, genesJ) ))
    df[i,j] = 1-jaccards
  }
}

## Cluster nodes using dynamic tree cut, for colors
distMat = as.dist(df)
dendro = hclust(distMat, method = "average" )
clust = dynamicTreeCut::cutreeDynamicTree( dendro, minModuleSize = 4 )


gs_results$Cluster = clust


df = df[,colnames(df) %in% gs_results$pathway]
df = df[rownames(df) %in% gs_results$pathway,]

distMat = as.dist(df)
dendro = hclust(distMat, method = "average" )

tal <- treeAndLeaf(dendro)

tal <- att.mapv(g = tal, dat = gs_results, refcol = 1)

tal <- att.setv(g = tal, from = "Cluster", to = "nodeColor", 
                cols = c('blue','red'))
tal <- att.setv(g = tal, from = "size", to = "nodeSize")

#--- Set graph attributes using 'att.addv' and 'att.adde' functions
tal <- att.addv(tal, "nodeFontSize", value = 15, index = V(tal)$isLeaf)
tal <- att.adde(tal, "edgeWidth", value = 3)

rdp <- RedPort()
calld(rdp)
resetd(rdp)
#--- Send the tree-and-leaf to the interactive R/Java interface
addGraph(obj = rdp, g = tal, gzoom=50)



addLegend.color(obj = rdp, tal, title = "Group", 
                position = "topright",vertical=T)

