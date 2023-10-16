### Script to separate get features of linnarson by Region to further classification with machine learning 

library('Matrix')
library("Seurat")
library('dplyr')
library("stringr")
library('ggplot2')
library("DUBStepR")

### Import of spared matrix, meta and genes of linnarson 

linnarson <- Matrix::readMM("../data/linnarson_sparse.mtx")
linnarson <- t(linnarson)

linnarson_meta <- read.csv("../data/linnarson_obs.csv",
                           header = T, stringsAsFactors = F)

linnarson_meta <- linnarson_meta[,c('CellID','Age','ClusterName', 'Clusters','Region', 'SampleID',
                                    'Sex', 'batch','MitoRiboRatio')]

linnarson_meta_var <- read.csv("../data/linnarson_var.csv",
                               header = T, stringsAsFactors = F)

### Get orthologus

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

### Creaete Seurat object

linnarson_seurat <- Seurat::CreateSeuratObject(linnarson_sub)

linnarson_seurat@meta.data$'CellID' <- colnames(linnarson_sub)
linnarson_seurat@meta.data$'Age' <- linnarson_meta$Age
linnarson_seurat@meta.data$'Region' <- linnarson_meta$Region
linnarson_seurat@meta.data$'Sex' <- linnarson_meta$Sex
linnarson_seurat@meta.data$'Cell_Type' <- linnarson_meta$batch
linnarson_seurat[["percent.mt"]] <- linnarson_meta$MitoRiboRatio


metadata <- linnarson_seurat@meta.data

# Rename columns
metadata <- metadata %>%
  dplyr::rename(
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata %>% 
  ggplot(aes(x=Cell_Type, fill=Cell_Type)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

metadata %>% 
  ggplot(aes(color=Cell_Type, x=nUMI, fill= Cell_Type)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=Cell_Type, x=nGene, fill= Cell_Type)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

metadata %>% 
  ggplot(aes(color=Cell_Type, x=percent.mt, fill=Cell_Type)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

linnarson_seurat$'log10GenesPerUMI' <- log10(linnarson_seurat$nUMI) / log10(linnarson_seurat$nGene)

### Filtered and subset

filtered_seurat <- subset(x = linnarson_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (percent.mt < 0.50) &
                             (Cell_Type != 'astrocytes') &
                            (Cell_Type != 'enteric_neurons'))

metadata <- filtered_seurat@meta.data

### rename areas with empty region in original metadata
metadata <- metadata %>% mutate(Region = case_when(
  Region == "" & Cell_Type == 'dentate_gyrus_granule_neurons' ~ 'Hippocampus',
  Region == "" & Cell_Type == 'olfactory_inhibitory_neurons' ~ 'Cortex',
  Region == "" & Cell_Type == 'peripheral_sensory_non_peptidergic_neurons' ~ 'PNS',
  Region == "" & Cell_Type == 'peripheral_sensory_peptidergic_neurons' ~ 'PNS',
  Region == "" & Cell_Type == 'sympathetic_noradrenergic_neurons' ~ 'PNS',
  TRUE ~ Region
))



  
filtered_seurat@meta.data <- metadata

### Subset cells based on excitatory, inhibitory and neurons from peripheric neuronal system
filtered_seurat <- subset(x = filtered_seurat,
                          subset = (Cell_Type %in% c('di_and_mesencephalon_excitatory_neurons','di_and_mesencephalon_inhibitory_neurons',
                                  'olfactory_inhibitory_neurons','peripheral_sensory_non_peptidergic_neurons','peripheral_sensory_peptidergic_neurons',
                                  'spinal_cord_excitatory_neurons','spinal_cord_inhibitory_neurons','sympathetic_noradrenergic_neurons',
                                  'telencephalon_inhibitory_interneurons','telencephalon_projecting_excitatory_neurons',
                                  'telencephalon_projecting_inhibitory_neurons')))

### create a column with excitatory, inhibitory and Peripherical neurons
filtered_seurat@meta.data <- filtered_seurat@meta.data %>% mutate(Neuron_Type = case_when(
  str_detect(Cell_Type,'excitatory') == TRUE ~ 'EXC',
  str_detect(Cell_Type,'inhibitory') == TRUE ~ 'INH',
  Cell_Type == 'di_and_mesencephalon_inhibitory_neurons' ~ 'INH',
  TRUE ~ 'PN'
))

### Change the name of the region to faclita visualization 

filtered_seurat@meta.data <- filtered_seurat@meta.data %>% mutate(Region =case_when(
  Region == 'Amygdala' ~ 'AMY',
  Region == 'Cortex' ~ 'CTX',
  Region %in% c('Dentate gyrus','Hippocampus','Hippocampus,Cortex') == TRUE ~ 'HIP',
  Region == 'Hypothalamus' ~ 'HYP',
  Region == 'Medulla' ~ 'MED',
  Region == 'Midbrain dorsal' ~ 'MDD',
  Region %in% c('Midbrain dorsal,Midbrain ventral','Midbrain ventral') == TRUE ~ 'MDV',
  Region == 'Pons' ~ 'PON' ,
  Region == 'Spinal cord' ~ 'SPC',
  Region == 'Striatum dorsal' ~ 'STD',
  Region == 'Striatum ventral' ~ 'STV',
  Region == 'Thalamus' ~ 'THA',
  TRUE ~ 'PNS'
))


### Remove cells from aged mice 

filtered_seurat <- subset(filtered_seurat,Age != 'p60')

### Get excitatory and peripherical neurons  and run DUBStepR

Region_EXC_PNS <- subset(filtered_seurat,Neuron_Type %in% c('EXC','PN'))

region_EXC <- data.frame(table(Region_EXC_PNS$Region))
region_EXC <- region_EXC %>% filter(Freq > 250) ### get only regions with more than 250 cells 


Region_EXC_PNS <- subset(Region_EXC_PNS, Region %in% region_EXC$Var1)

Region_EXC_PNS_RNA <- NormalizeData(object = Region_EXC_PNS, normalization.method = "LogNormalize")
Idents(Region_EXC_PNS_RNA) <- 'Region'

dubstepR.out_EXC <- DUBStepR(input.data = Region_EXC_PNS_RNA@assays$RNA@data,
                             min.cells = 0.05*ncol(Region_EXC_PNS_RNA), 
                             optimise.features = T, k = 10, num.pcs = 20, error = 0)
 

all_data_DUB_EXC <- readRDS("../data/all_data_EXC_DUB.rds")
Idents(all_data_DUB_EXC) <- 'RNA'

all_data_DUB_EXC <- as.data.frame(GetAssayData(all_data_DUB_EXC))[,]
all_data_DUB_EXC <- as.data.frame(t(all_data_DUB_EXC))
all_data_DUB_EXC <- all_data_DUB_EXC[,order(colnames(all_data_DUB_EXC))]
all_data_DUB_EXC <- all_data_DUB_EXC[,colnames(all_data_DUB_EXC) %in% dubstepR.out_EXC$optimal.feature.genes]
write.csv(all_data_DUB_EXC,'../data/all_data_Region_EXC.csv')

DefaultAssay(Region_EXC_PNS) <- 'RNA'
Region_EXC_PNS_region <- Region_EXC_PNS$Region

Region_EXC_PNS <- as.data.frame(GetAssayData(Region_EXC_PNS))[,]
Region_EXC_PNS <- as.data.frame(t(Region_EXC_PNS))
Region_EXC_PNS <- Region_EXC_PNS[,order(colnames(Region_EXC_PNS))]
Region_EXC_PNS <- Region_EXC_PNS[,colnames(Region_EXC_PNS) %in% colnames(all_data_DUB_EXC)]
Region_EXC_PNS$'Region' <- Region_EXC_PNS_region
write.csv(Region_EXC_PNS,'../data/Region_EXC_PNS.csv')



### Get inhibitory and peripherical neurons  and run DUBStepR

Region_INH_PNS <- subset(filtered_seurat,Neuron_Type %in% c('INH','PN'))

region_INH <- data.frame(table(Region_INH_PNS$Region))
region_INH <- region_INH %>% filter(Freq > 250) ### get only regions with more than 250 cells 


Region_INH_PNS <- subset(Region_INH_PNS, Region %in% region_INH$Var1)

Region_INH_PNS_RNA <- NormalizeData(object = Region_INH_PNS, normalization.method = "LogNormalize")
Idents(Region_INH_PNS_RNA) <- 'Region'

dubstepR.out_INH <- DUBStepR(input.data = Region_INH_PNS_RNA@assays$RNA@data,
                             min.cells = 0.05*ncol(Region_INH_PNS_RNA), 
                             optimise.features = T, k = 10, num.pcs = 20, error = 0)

all_data_DUB_INH <- readRDS("../data/all_data_INH_DUB.rds")
Idents(all_data_DUB_INH) <- 'RNA'

all_data_DUB_INH <- as.data.frame(GetAssayData(all_data_DUB_INH))[,]
all_data_DUB_INH <- as.data.frame(t(all_data_DUB_INH))
all_data_DUB_INH <- all_data_DUB_INH[,order(colnames(all_data_DUB_INH))]
all_data_DUB_INH <- all_data_DUB_INH[,colnames(all_data_DUB_INH) %in% dubstepR.out_INH$optimal.feature.genes]
write.csv(all_data_DUB_INH,'../data/all_data_Region_INH.csv')

DefaultAssay(Region_INH_PNS) <- 'RNA'
Region_INH_PNS_region <- Region_INH_PNS$Region

Region_INH_PNS <- as.data.frame(GetAssayData(Region_INH_PNS))[,]
Region_INH_PNS <- as.data.frame(t(Region_INH_PNS))
Region_INH_PNS <- Region_INH_PNS[,order(colnames(Region_INH_PNS))]
Region_INH_PNS <- Region_INH_PNS[,colnames(Region_INH_PNS) %in% colnames(all_data_DUB_INH)]
Region_INH_PNS$'Region' <- Region_INH_PNS_region
write.csv(Region_INH_PNS,'../data/Region_INH_PNS.csv')
