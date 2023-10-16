library("Seurat")
library("stringr")
library("Matrix")
library("clustifyr")
library("dplyr")
### Integration of SEA AD dataset (https://cellxgene.cziscience.com/collections/1ca90a2d-2943-483d-b678-b809bf464c30)

### import of metadata

SEAD_meta <- read.csv("SEA_AD/SEAD_metadata.csv",header = T, stringsAsFactors = F)
SEAD_meta <- SEAD_meta[SEAD_meta$Diagnosis %in% c('Alzheimers disease','Control'),]


### Subset the seurat objects

seurat_objects<- list.files("SEA_AD/SEAD_objects/")

seurat_objects <- paste0("SEA_AD/SEAD_objects/",seurat_objects)


seurat_names <- as.data.frame(str_split_fixed(seurat_objects,"/",3))
seurat_names <- seurat_names$V3
seurat_names <- str_remove(seurat_names,".rds")


for (i in 1:21) {

  seurat <- readRDS(seurat_objects[i])
  seurat@meta.data$'cell_id' <- colnames(seurat)

  meta_seurat <- seurat@meta.data
  meta_seurat <- meta_seurat[meta_seurat$donor_id %in% SEAD_meta$donor_id,]
  meta_seurat <- merge(meta_seurat,SEAD_meta,by='donor_id')
  meta_seurat <- meta_seurat %>% group_by(Diagnosis) %>% slice_sample(n=1250) ### Minimal cell by dataset

  seurat@meta.data$'cells_to_keep' <- seurat$cell_id %in% meta_seurat$cell_id

  seurat <- subset(seurat,cells_to_keep == TRUE)
  print(seurat_objects[i])


  saveRDS(seurat,paste0("SEAD_objects/",seurat_names[i],".rds"))

  rm(seurat)

}


scrna.list <- list.files("SEA_AD/SEAD_objects/")
scrna.list <- paste("SEAD_objects",scrna.list,sep = '/')

scrna_objects_list <- list()

for (i in 1:21) {

  seurat <- readRDS(scrna.list[i])
  scrna_objects_list[[i]] <- seurat
  rm(seurat)

}

scrna <- merge(x=scrna_objects_list[[1]], y=c(scrna_objects_list[[2]],scrna_objects_list[[3]],
                                       scrna_objects_list[[4]],
                                       scrna_objects_list[[5]],
                                       scrna_objects_list[[6]],
                                       scrna_objects_list[[7]],
                                       scrna_objects_list[[8]],
                                       scrna_objects_list[[9]],
                                       scrna_objects_list[[10]],
                                       scrna_objects_list[[11]],
                                     scrna_objects_list[[12]],
                                      scrna_objects_list[[13]],
                                      scrna_objects_list[[14]],
                                      scrna_objects_list[[15]],
                                      scrna_objects_list[[16]],
                                       scrna_objects_list[[17]],
                                       scrna_objects_list[[18]],
                                       scrna_objects_list[[19]],
                                      scrna_objects_list[[20]],
                                     scrna_objects_list[[21]]),project="SEAD",
 add.cell.ids = c("A","B","C",'D','E','F','G','H','I','J','L','M','N','O','P','Q','R','S','T','U','V'))


features_name <- data.frame(genes = scrna_objects_list[[1]]@assays$RNA@meta.features$feature_name,stringsAsFactors = F)
features_name$genes <- as.character(features_name$genes)

scrna@assays$RNA@counts@Dimnames[[1]] <- features_name$genes
scrna@assays$RNA@data@Dimnames[[1]] <- features_name$genes
scrna@assays$RNA@meta.features <- features_name

scrna_meta <- scrna@meta.data

SEAD_meta <- read.csv("SEAD_metadata.csv",header = T, stringsAsFactors = F)
SEAD_meta <- SEAD_meta[SEAD_meta$Diagnosis %in% c('Alzheimers disease','Control'),]
SEAD_meta <- SEAD_meta[SEAD_meta$donor_id %in% scrna_meta$donor_id,]
SEAD_meta <- merge(SEAD_meta,scrna_meta,by='donor_id')
SEAD_meta <- SEAD_meta[order(match(SEAD_meta$cell_id,scrna_meta$cell_id)),]

scrna@meta.data$'Diagnosis' <- SEAD_meta$Diagnosis
saveRDS(scrna,"SEAD_sub.rds")

