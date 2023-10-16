###3 Feature selection of linnarson with DUBStepR

library("Seurat")
library("scrattch.io")
library("edgeR")
library("dplyr")
library("DUBStepR")
library("lognorm")
library("scCustomize")
library('ggplot2')



linnarson <- readRDS("~/Análise/linnarson/linnarson_with_astrocytes/linnarson_seurat_sub.rds")


### Percentage of HBN (Habenullar Neurons) that are GABAergic

linnarson[['GAD1']] <- PercentageFeatureSet(linnarson,pattern = "GAD1")


### Classification of cell types in glutamatergic or GABAergic 

linnarson@meta.data <- linnarson@meta.data %>% mutate(cell_type = case_when(
Cell == 'HBN' & GAD1 > 0.0 ~ 'INH',
Cell == 'HBN' & GAD1 == 0.0 ~'EXC',
Cell %in% c('DGGN','DMEN','SCEN','TPEN','CBN') == T ~ 'EXC',
Cell %in% c('DMIN','OIN','SCIN','TII','TPIN') == T ~ 'INH',
TRUE ~ 'exclude!'
))


### Exclusion of unwanted cells 

linnarson <- subset(linnarson_meta,cell_type != 'exclude!')

### Basice pipeline of DUBStepR

DefaultAssay(linnarson) <- 'RNA'

linnarson <-  NormalizeData(object = linnarson, normalization.method = "LogNormalize")

dubstepR.out <- DUBStepR(input.data = linnarson@assays$RNA@data, min.cells = 0.05*ncol(linnarson), 
                         optimise.features = T, k = 10, num.pcs = 20, error = 0)
linnarson@assays$RNA@var.features <- dubstepR.out$optimal.feature.genes
linnarson <- ScaleData(linnarson, features = rownames(linnarson))
linnarson_meta <- RunPCA(linnarson, features =  VariableFeatures(object = linnarson_meta))


seuratObj <- FindNeighbors(linnarson, reduction = "pca", dims = 1:50)
seuratObj <- FindClusters(seuratObj, resolution = c(0.1,0.2,0.4,0.8,1,1.2))
seuratObj <- RunUMAP(seuratObj, dims = 1:50, n.components = 2)


features_corr <- dubstepR.out$optimal.feature.genes

Idents(seuratObj) <- 'cell_type'

### Find top markers based on cell type

top.markers <- FindAllMarkers(object = seuratObj, assay = "RNA", logfc.threshold = 0.25, min.pct = 0.2, only.pos = FALSE) 
top.markers <- top.markers %>% filter(p_val_adj < 0.01)

### Save objects

saveRDS(top.markers,'Análise/linnarson/linnarson_with_astrocytes/top_markers_DUB.rds')
saveRDS(seuratObj,'Análise/linnarson/linnarson_with_astrocytes/linnarson_seurat_DUB.rds')



### Write subset of matrix with top markers 
DefaultAssay(linnarson_meta) <- 'RNA'

mat_all <- GetAssay(linnarson_meta)
mat_all <- mat_all@counts
mat_all <- mat_all[top.markers$gene,]
write_dgCMatrix_csv(mat_all,"Análise/linnarson/linnarson_with_astrocytes/linnarson_sub_DUB.csv",col1_name = 'gene')













