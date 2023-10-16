### Import of seurat objects - Only Neurons

### Integration of all datasets 


library("Seurat")
library("dplyr")
library('RCurl')
library('harmony')
library('ggplot2')
## ASCL1

ASCL1 <- readRDS("only_neurons/Ascl1_iNs.rds")
DefaultAssay(ASCL1) <- 'RNA'
ASCL1 <- CreateSeuratObject(counts = ASCL1@assays$RNA@counts)
ASCL1$'dataset' <- 'ASCL1'



### hiNPC - 4 and 6 weeks

hiNPCs_4w <- readRDS("only_neurons/hiNPCs_4w.rds")
DefaultAssay(hiNPCs_4w) <- 'RNA'
hiNPCs_4w <- CreateSeuratObject(counts = hiNPCs_4w@assays$RNA@counts)
hiNPCs_4w$'dataset' <- 'hiNPCs_4w'

hiNPCs_6w <- readRDS("only_neurons/hiNPCs_6w.rds")
DefaultAssay(hiNPCs_6w) <- 'RNA'
hiNPCs_6w <- CreateSeuratObject(counts = hiNPCs_6w@assays$RNA@counts)
hiNPCs_6w$'dataset' <- 'hiNPCs_6w'

### NEUROG2

NEUROG2 <- readRDS("only_neurons/Neurog2_iNs.rds")

DefaultAssay(NEUROG2) <- 'RNA'
NEUROG2 <- CreateSeuratObject(counts = NEUROG2@assays$RNA@counts)
NEUROG2$'dataset' <- 'NEUROG2'


#### Organoids

Organoids <- readRDS("only_neurons/org_iNs.rds")

DefaultAssay(Organoids) <- 'RNA'
Organoids <- CreateSeuratObject(counts = Organoids@assays$RNA@counts)
Organoids$'dataset' <- 'Organoids'


seurat_list <- list(ASCL1,NEUROG2,hiNPCs_4w,hiNPCs_6w,Organoids)

# Merge raw samples

merged_seurat <- merge(x = seurat_list[[1]],
                       y = seurat_list[2:length(seurat_list)],
                       merge.data = TRUE)

merged_seurat@meta.data <- merged_seurat@meta.data %>% select(orig.ident,nCount_RNA,nFeature_RNA,dataset)
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")

### Regress cell cycle genes
library('AnnotationHub')

cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") 
cell_cycle_genes <- read.csv(text = cc_file)

ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

### Calculate th cell scoring genes

merged_seurat <- CellCycleScoring(merged_seurat,
                                  g2m.features = g2m_genes,
                                  s.features = s_genes)

### Normalize and run SCT with percent.mt, S.Score and G2M.Scores as covarietes

merged_seurat <- merged_seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>% 
  ScaleData() %>%
  SCTransform(vars.to.regress = c("percent.mt",'S.Score','G2M.Score'))



merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 30)


### Integrate merged object with Harmony 

harmonized_seurat <- RunHarmony(merged_seurat, 
                                group.by.vars = c("dataset"), 
                                reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:30)

harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = c(0.1,0.2, 0.4, 0.6, 0.8, 1.0, 1.2))
Idents(harmonized_seurat) <- "SCT_snn_res.0.4"

### Save
saveRDS(harmonized_seurat,"only_neurons_tables/harmonized_all_seurat.rds")


