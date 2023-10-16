library("Seurat")
library("scrattch.io")
library("edgeR")
library("dplyr")
library("DUBStepR")
library("lognorm")
library("scCustomize")
library('ggplot2')
library('dplyr')

### 

SEAD_sub <- readRDS('~/SEA_AD/SEAD_sub.rds')

class_to_keep <- c('Chandelier','L4 IT','L5 ET','L5 IT','L5/6 NP','L6 CT','L6 IT',
                  'L6 IT Car3','L6b','Lamp5','Lamp5 Lhx6','Pax6','Pvalb','Sncg','Sst','Vip') ### Select only excitatory and inhibitory neurons

SEAD_sub <- subset(SEAD_sub, Subclass %in% class_to_keep )
SEAD_sub <- subset(SEAD_sub, Diagnosis == 'Control') ### get only cells from control patients
DefaultAssay(SEAD_sub) <- 'RNA'

#### Subset SEAD eith the same genes used in classification of excitatory and inhibitory neurons with linnarson dataset

DUB_markers <- read.csv("../data/linnarson_sub_DUB.csv",header = T, stringsAsFactors = F,row.names = 1)
DUB_markers <- colnames(DUB_markers)
DUB_markers <- DUB_markers[1:533]

SEAD_sub <- as.data.frame(GetAssayData(SEAD_sub))[,]
SEAD_sub <- as.data.frame(t(SEAD_sub))
SEAD_sub <- SEAD_sub[,order(colnames(SEAD_sub))]
SEAD_sub <- all_data_DUB_EXC[,colnames(SEAD_sub) %in% DUB_markers]

write.csv(SEAD_sub,"../data/SEAD_sub.csv")





