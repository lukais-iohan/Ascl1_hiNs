### Script to make fluvial plots with prediction of regions 

library('ggalluvial')
library("dplyr")
library('ggsankey')
library('Seurat')

### Import of necessary data

all_data_EXC <- readRDS("../data/all_data_EXC_DUB.rds")
all_data_INH <- readRDS("../data/all_data_INH_DUB.rds")

predicted_EXC_Region <- read.csv("../data/all_data_sensory_EXC_predicted.csv",header = T, 
                                 stringsAsFactors = F, row.names = 1)

predicted_INH_Region <- read.csv("../data/all_data_sensory_INH_.csv",header = T, 
                                 stringsAsFactors = F, row.names = 1)

all_data_EXC <- all_data_EXC@meta.data
all_data_EXC$'SVM' <- predicted_EXC_Region$SVM
all_data_EXC$'MLP' <- predicted_EXC_Region$MLP

all_data_INH <- all_data_INH@meta.data
all_data_INH$'SVM' <- predicted_INH_Region$SVM
all_data_INH$'MLP' <- predicted_INH_Region$MLP

### Function to plot fluvial EXC
fluvial_function_EXC <- function(table,AI_algorithm,neuron_type) {

  all_data_meta <- data.frame(dataset = table$dataset,
                              AI_Method = rep(AI_algorithm,nrow(table)), 
                              Neuron = table[,AI_algorithm])
  
  TotalCount = nrow(all_data_meta)
  # Step 1
  all_data_longer <- all_data_meta %>%
    make_long( dataset, Neuron)
  
  # Step 2
  dagg <- all_data_longer%>%
    dplyr::group_by(node)%>%
    tally()
  
  dagg <- dagg%>%
    dplyr::group_by(node)%>%
    dplyr::mutate(pct = n/TotalCount)
  
  
  # Step 3
  df2 <- merge(all_data_longer, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)
  
  df2$node <- factor(df2$node, 
  levels = c('AMY','ASCL1','CTX','NEUROG2','hiNPCs_6w','HIP','MDD','MED','hiNPCs_4w','PNS','SPC','THA'))
  
  pl <- ggplot(df2, aes(x = x
                        , next_x = next_x
                        , node = node
                        , next_node = next_node
                        , fill = factor(node)
                        
                        , label = paste0(node," n=", n, '(',  round(pct* 100,1), '%)' ))
  )
  
  pl + geom_sankey(flow.alpha = 0.6, node.color = 1) +
    geom_sankey_label(size = 5, color = 1, fill = "white") +
    theme_sankey(base_size = 16) +scale_fill_hue()+
    theme(legend.position = "none") + ggtitle(label = paste0('Prediction with ', AI_algorithm,' ', neuron_type)) +xlab('')
   
  
}

### EXC fluvial plot

lapply('SVM', function(i) {
  
  fluvial_function_EXC(all_data_EXC,i,'EXC')
})

### Figure saved with rstudio interface

### tables of EXC frequency of each region

library('gt')

MLP_EXC <- data.frame(table(all_data_EXC$dataset,all_data_EXC$SVM))

MLP_EXC <- MLP_EXC %>% group_by(Var1) %>%
  mutate(freq = (round(Freq / sum(Freq)*100)))

MLP_EXC <- MLP_EXC %>% select(Var1,Var2,freq)

colnames(MLP_EXC) <- c('Dataset','Region','Frequency')


MLP_EXC  %>% group_by(Dataset) %>%
    gt()
  

### Function to plot INH fluvial

fluvial_function_INH <- function(table,AI_algorithm,neuron_type) {
  
  all_data_meta <- data.frame(dataset = table$dataset,
                              AI_Method = rep(AI_algorithm,nrow(table)), 
                              Neuron = table[,AI_algorithm])
  
  TotalCount = nrow(all_data_meta)
  # Step 1
  all_data_longer <- all_data_meta %>%
    make_long( dataset, Neuron)
  
  # Step 2
  dagg <- all_data_longer%>%
    dplyr::group_by(node)%>%
    tally()
  
  dagg <- dagg%>%
    dplyr::group_by(node)%>%
    dplyr::mutate(pct = n/TotalCount)
  
  
  # Step 3
  df2 <- merge(all_data_longer, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)
  
  df2$node <- factor(df2$node, 
                     levels = c('ASCL1','CTX','NEUROG2','hiNPCs_6w','HIP','MDD','MDV','MED','hiNPCs_4w','PNS',
                                'SPC','STD','STV','THA'))
  
  pl <- ggplot(df2, aes(x = x
                        , next_x = next_x
                        , node = node
                        , next_node = next_node
                        , fill = factor(node)
                        
                        , label = paste0(node," n=", n, '(',  round(pct* 100,1), '%)' ))
  )
  
  pl + geom_sankey(flow.alpha = 0.6, node.color = 1) +
    geom_sankey_label(size = 5, color = 1, fill = "white") +
    theme_sankey(base_size = 16) +scale_fill_hue()+
    theme(legend.position = "none") + ggtitle(label = paste0('Prediction with ', AI_algorithm,' ', neuron_type)) +xlab('')
  
  
}

lapply('SVM', function(i) {
  
  fluvial_function_INH(all_data_INH,i,'INH')
})

### Figure saved with rstudio interface

### table of INH frequency of each region

MLP_INH <- data.frame(table(all_data_INH$dataset,all_data_INH$SVM))

MLP_INH <- MLP_INH %>% group_by(Var1) %>%
  mutate(freq = (round(Freq / sum(Freq)*100)))

MLP_INH <- MLP_INH %>% select(Var1,Var2,freq)

colnames(MLP_INH) <- c('Dataset','Region','Frequency')


MLP_INH  %>% group_by(Dataset) %>%
  gt()

