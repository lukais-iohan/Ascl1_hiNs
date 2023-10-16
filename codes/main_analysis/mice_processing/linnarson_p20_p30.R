### get only animals from p20-p30 in linnarson sparse

library('Matrix')
library("Seurat")
library('dplyr')
library("stringr")

### Import of sparsed matrix, metadata and genes from linnarson

linnarson <- Matrix::readMM("~/Análise/linnarson/linnarson_with_astrocytes/linnarson_sparse.mtx")
linnarson <- t(linnarson)

linnarson_meta <- read.csv("~/Análise/linnarson/linnarson_with_astrocytes/linnarson_obs.csv",
                            header = T, stringsAsFactors = F)

linnarson_meta <- linnarson_meta[,c('CellID','Age','ClusterName', 'Clusters','Region', 'SampleID',
                                     'Sex', 'batch','MitoRiboRatio')]

linnarson_meta_var <- read.csv("~/Análise/linnarson/linnarson_with_astrocytes/linnarson_var.csv",
                                header = T, stringsAsFactors = F)

### Get ortohlogus of human genes with Gprofiler

Ort <- gprofiler2::gorth(linnarson_meta_var$Gene,source_organism = 'mmusculus',target_organism =  'hsapiens')
Ort <- na.omit(Ort)

linnarson@Dimnames[[1]] <- linnarson_meta_var$Gene
linnarson@Dimnames[[2]] <- linnarson_meta$CellID

linnarson_sub <- linnarson[which(rownames(linnarson)%in% Ort$input),]

Ort <- Ort[rownames(linnarson_sub) %in% Ort$input,]
Ort <- Ort[!duplicated(Ort$input),]
Ort_sorted <- Ort[order(match(x=rownames(linnarson_sub),table = Ort$input)),]

linnarson_sub@Dimnames[[1]] <- Ort$ortholog_name
linnarson_sub@Dim <- as.integer(c('16960','88234'))


### Create Seurat Object

linnarson_seurat <- Seurat::CreateSeuratObject(linnarson_sub)

linnarson_seurat@meta.data$'CellID' <- colnames(linnarson_sub)
linnarson_seurat@meta.data$'Age' <- linnarson_meta$Age
linnarson_seurat@meta.data$'Region' <- linnarson_meta$Region
linnarson_seurat@meta.data$'Sex' <- linnarson_meta$Sex
linnarson_seurat@meta.data$'Cell_Type' <- linnarson_meta$batch

linnarson_seurat[["percent.mt"]] <- linnarson_meta$MitoRiboRatio


### Filtering of low quality cells

linnarson_seurat_sub <- subset(linnarson_seurat, subset = nFeature_RNA > 200
                            & nFeature_RNA < 5000 & percent.mt < 5 & nCount_RNA < 20000)


### Change the names of cell_type

linnarson_seurat_sub@meta.data <- linnarson_seurat_sub@meta.data %>% mutate(Cell = case_when(
  linnarson_seurat_sub@meta.data$Cell_Type == 'astrocytes' ~ 'AST',
  linnarson_seurat_sub@meta.data$Cell_Type == 'cerebellum_neurons' ~ 'CBN',
  linnarson_seurat_sub@meta.data$Cell_Type == 'cholinergic_and_monoaminergic_neurons' ~ 'CMN',
  linnarson_seurat_sub@meta.data$Cell_Type == 'dentate_gyrus_granule_neurons' ~ 'DGGN',
  linnarson_seurat_sub@meta.data$Cell_Type == 'di_and_mesencephalon_excitatory_neurons' ~ 'DMEN',
  linnarson_seurat_sub@meta.data$Cell_Type == 'di_and_mesencephalon_inhibitory_neurons' ~ 'DMIN',
  linnarson_seurat_sub@meta.data$Cell_Type == 'enteric_neurons' ~ 'EN',
  linnarson_seurat_sub@meta.data$Cell_Type == 'hindbrain_neurons' ~ 'HBN',
  linnarson_seurat_sub@meta.data$Cell_Type == 'olfactory_inhibitory_neurons' ~ 'OIN',
  linnarson_seurat_sub@meta.data$Cell_Type == 'peptidergic_neurons' ~ 'PPN',
  linnarson_seurat_sub@meta.data$Cell_Type == 'peripheral_sensory_non_peptidergic_neurons' ~ 'PSNPN',
  linnarson_seurat_sub@meta.data$Cell_Type == 'peripheral_sensory_peptidergic_neurons' ~ 'PSPN',
  linnarson_seurat_sub@meta.data$Cell_Type == 'spinal_cord_excitatory_neurons' ~ 'SCEN',
  linnarson_seurat_sub@meta.data$Cell_Type == 'spinal_cord_inhibitory_neurons' ~ 'SCIN',
  linnarson_seurat_sub@meta.data$Cell_Type == 'sympathetic_noradrenergic_neurons' ~ 'SNAN',
  linnarson_seurat_sub@meta.data$Cell_Type == 'telencephalon_inhibitory_interneurons' ~ 'TII',
  linnarson_seurat_sub@meta.data$Cell_Type == 'telencephalon_projecting_excitatory_neurons' ~ 'TPEN',
  TRUE ~ 'TPIN'

))




### Subset to get only animals p20-p30
linnarson_seurat_sub <- subset(linnarson_seurat_sub,Age %in% c('p20','p21','p22','p25-27','p25','p26','p27','p28','p29','p30'))

### Normalizing by percent.mt and save

 linnarson_seurat_sub <- linnarson_seurat_sub %>% 
  SCTransform(vars.to.regress = "percent.mt") %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

saveRDS(linnarson_seurat_sub,"~/Análise/linnarson/linnarson_with_astrocytes/linnarson_seurat_sub.rds")
