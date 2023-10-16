library("Seurat")


#### Import of EXC and INH seurat objects

all_data_EXC_DUB <- readRDS("..data/all_data_EXC_DUB.rds")
all_data_INH_DUB <- readRDS("../data/all_data_INH_DUB.rds")

Idents(all_data_EXC_DUB) <- 'dataset'
Idents(all_data_INH_DUB) <- 'dataset'

### Function to find differential expression genes

DE_Seurat <- function(seurat_object,dataset1,dataset2) {
  
  
  markers <- FindMarkers(seurat_object,ident.1 = dataset1,ident.2 = dataset2,min.pct = 0.2,
                         logfc.threshold = 0.25)
  
  markers$'comparasion' <- paste(dataset1,dataset2,sep= " x ")
  markers$"genes" <- rownames(markers)
  return(markers)
  
}


### Function to compare make all possible combinations between datasets

all_data_DEGs <- function(dataset) {
  
  data1 <- c('ASCL1','hiNPCs_6w','NEUROG2','hiNPCs_4w')
  data2 <- c('ASCL1','hiNPCs_6w','NEUROG2','hiNPCs_4w')
  
  seurat_degs <- list()
  
  for (i in 1:(length(data1)-1)) {
    for (j in (i+1):length(data1)) {
      seurat_degs[[paste(data1[i], data2[j], sep = "_")]] <- DE_Seurat(dataset, dataset1 = data1[i], 
                                                                       dataset2 = data2[j])
    }
  }
  return(seurat_degs)
}


### Get all comparasion tables to only one table to both EXC and INH

EXC_DEGs <- all_data_DEGs(all_data_EXC_DUB)
EXC_DEGs <- plyr::join_all(EXC_DEGs,type = 'full')
EXC_DEGs$'metric' <-  sign(EXC_DEGs$avg_log2FC)*rank(abs(EXC_DEGs$avg_log2FC)*-log10(EXC_DEGs$p_val_adj)) ### metric to run GSEA 

INH_DEGs <- all_data_DEGs(all_data_INH_DUB)
INH_DEGs <- plyr::join_all(INH_DEGs,type = 'full')
INH_DEGs$'metric' <-  sign(INH_DEGs$avg_log2FC)*rank(abs(INH_DEGs$avg_log2FC)*-log10(INH_DEGs$p_val_adj))  ### metric to run GSEA 

write.csv(EXC_DEGs,'only_neurons_tables/EXC_DEGs.csv')
write.csv(INH_DEGs,'only_neurons_tables/INH_DEGs.csv')
