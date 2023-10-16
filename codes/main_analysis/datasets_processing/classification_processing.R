
#### Classification processing

library("Seurat")


### Import of seurat object with ASCL1, hiNPCs_4w, hiNPCs_6w, NEUROG2 and Organoids integrated data

all_data <- readRDS("../data/harmonized_all_seurat.rds")


### Import of prediciton cell type results
all_data_predicted <- read.csv("../data/all_data_DUB_grid_predicted.csv",
                               header = T, stringsAsFactors = F, row.names = 1)



Idents(all_data) <- 'dataset'


### Get only cells with the same classification by SVM and MLP algortihms and exclude Organoids cells (not include in the analysis)
all_data$'MLP' <- all_data_predicted$MLP

all_data$'SVM' <- all_data_predicted$SVM

all_data <- subset(all_data,subset = dataset != 'Organoids')

all_data$'same_classification' <- all_data$SVM == all_data$MLP

all_data <- subset(all_data, same_classification != FALSE)

Idents(all_data) <- 'dataset_classification'


all_data$'dataset_classification' <- paste(all_data$dataset,all_data$SVM,sep='.')

saveRDS(all_data,'only_neurons_tables/all_data_DUB.rds')

##### Split by neuron type 

all_data_split <- SplitObject(all_data,split.by = 'MLP')


all_data_INH <- CreateSeuratObject(counts = all_data_split$INH@assays$RNA@counts)

all_data_INH$'dataset' <- all_data_split$INH$dataset

all_data_EXC <- CreateSeuratObject(counts = all_data_split$EXC@assays$RNA@counts)

all_data_EXC$'dataset' <- all_data_split$EXC$dataset

### Split datasets based on cell type classification and run DUBStepR

### INH


all_data_INH_RNA <- NormalizeData(object = all_data_INH, normalization.method = "LogNormalize")

dubstepR.out_INH <- DUBStepR(input.data = all_data_INH_RNA@assays$RNA@data,
                         min.cells = 0.05*ncol(all_data_INH), 
                         optimise.features = T, k = 10, num.pcs = 20, error = 0)
 

all_data_INH <- SCTransform(all_data_INH, vst.flavor = "v2", verbose = T,
                            residual.features = dubstepR.out_INH$optimal.feature.genes) 
  

all_data_INH <- RunPCA(all_data_INH, features =  VariableFeatures(object = all_data_INH))

all_data_INH <- FindNeighbors(all_data_INH,reduction = "pca", dims = 1:dubstepR.out_INH$elbow.pt, verbose = FALSE) 
all_data_INH <- FindClusters(all_data_INH,resolution = c(0.1,0.2,0.4,0.8,1,1.2), verbose = FALSE)
all_data_INH <- RunUMAP(all_data_INH, dims = 1:dubstepR.out_INH$elbow.pt, verbose = FALSE) 


DimPlot_scCustom(seurat_object = all_data_INH ,pt.size = 0.5,group.by = 'dataset',label = T,
                       ggplot_default_colors = T,figure_plot = T)
DimPlot_scCustom(seurat_object = all_data_INH ,pt.size = 0.5,group.by = 'SCT_snn_res.0.1',label = T,
                       ggplot_default_colors = T,figure_plot = T)


### EXC

all_data_EXC_RNA <- NormalizeData(object = all_data_EXC, normalization.method = "LogNormalize")

dubstepR.out_EXC <- DUBStepR(input.data = all_data_EXC_RNA@assays$RNA@data,
                         min.cells = 0.05*ncol(all_data_EXC), 
                         optimise.features = T, k = 10, num.pcs = 20, error = 0)

all_data_EXC <- SCTransform(all_data_EXC, vst.flavor = "v2", verbose = T,
                            residual.features = dubstepR.out_EXC$optimal.feature.genes) 

all_data_EXC<- RunPCA(all_data_EXC, features =  VariableFeatures(object = all_data_EXC))

all_data_EXC <- FindNeighbors(all_data_EXC,reduction = "pca", dims = 1:dubstepR.out_EXC$elbow.pt, verbose = FALSE) 
all_data_EXC <- FindClusters(all_data_EXC,resolution = c(0.1,0.2,0.4,0.8,1,1.2), verbose = FALSE)
all_data_EXC <- RunUMAP(all_data_EXC, dims = 1:dubstepR.out_EXC$elbow.pt, verbose = FALSE) 


### Save

saveRDS(all_data_EXC,'../data/all_data_EXC_DUB.rds')
saveRDS(all_data_INH,'../data/all_data_INH_DUB.rds')
