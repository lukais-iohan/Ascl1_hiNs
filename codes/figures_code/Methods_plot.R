#### Figure 1 script

library("Seurat")
library("plyr")
library('dplyr')
library('ggplot2')
library('hrbrthemes')
library('patchwork')
library('scCustomize')
library('grid')
library('viridis')
library('stringr')
library('cowplot')
library('RColorBrewer')
library('ggrepel')
library('dittoSeq')
library('cowplot')
library('viridis')
library('scales')
library('ggrepel')

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

dim_plot <- DimPlot_scCustom(seurat_object = all_data ,pt.size = 0.5,group.by = 'dataset',label = F,
                 ggplot_default_colors = T,repel = T)

dim_neuron <- DimPlot_scCustom(seurat_object = all_data ,pt.size = 0.5,group.by = 'SVM',label = F,
                 colors_use = c('firebrick2','dodgerblue2'),repel = T) 

fig1 <- dim_plot/dim_neuron

ggsave('figuras_artigo/fig1/dimplot.jpeg',plot = fig1,width = 8,height = 8,dpi = 600)

### Bar plot
dittoBarPlot(all_data,scale = 'percent',group.by = 'dataset',var = 'SVM',
             theme = theme_ipsum(base_size = 20,axis_text_size = 20,
                                 plot_title_size = 50,grid = F,axis_title_size = 20,axis_title_family = 20),
             color.panel = c('firebrick2','dodgerblue2')) + ggtitle('') + theme(axis.text = element_text(color = 'black'),
                        axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size = 25),
                        legend.text = element_blank()) + 
  xlab('') + ylab('')

ggsave('figuras_artigo/fig1/barplot.jpeg',width = 8,height = 6,units = 'in',dpi = 600)

#### Violin plot

gene_list_plot <- c("MAPT", "RBFOX3","SLC18A3","SLC6A3" , "SLC17A6", "GAD1", "GAD2")



Stacked_VlnPlot(seurat_object = all_data, features = gene_list_plot , x_lab_rotate = T,
                colors_use = c('firebrick2','dodgerblue2'), split.by = "SVM",plot_legend = T) 


ggsave('figuras_artigo/fig1/violin_plot.jpeg',width = 10,height = 5,dpi = 600)
dev.off()
 
