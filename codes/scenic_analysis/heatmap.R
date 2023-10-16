##### Heatmap SCENIC
library('tidyr')
library('dplyr')
library('plotrix')
library('ggplot2')
library('hrbrthemes')
library('stringr')

### Heatmap of regulon activity

EXC_SCENIC <- read.csv("results/regulon_EXC_dataset.csv",header = T, stringsAsFactors = F,row.names = 1)
EXC_SCENIC_scaled <-  t(scale((EXC_SCENIC), center = T, scale=T))

ComplexHeatmap::Heatmap(EXC_SCENIC_scaled,name = '-',cluster_rows = T,cluster_columns = F)


INH_SCENIC <- read.csv("results/regulon_INH_dataset_sampled.csv",header = T, stringsAsFactors = F,row.names = 1)
INH_SCENIC_scaled <-  t(scale((INH_SCENIC), center = T, scale=T))

ComplexHeatmap::Heatmap(INH_SCENIC_scaled,name = '-',cluster_rows = T,cluster_columns = F)



