### NEtwork of scenic regulons of EXC dataset

library('stringr')
library('visNetwork')
library('igraph')
library('RedeR')
library('scales')
library('dplyr')
#### EXC

adj_EXC <- read.csv('results/adj_EXC.csv',header = T, stringsAsFactors = F)


EXC_DEGs <- read.csv("../data/EXC.csv",header = T, stringsAsFactors = F)
EXC_DEGs <- EXC_DEGs[EXC_DEGs$diff_expressed != 'NO',]
EXC_DEGs <- EXC_DEGs[EXC_DEGs$comparasion %in% c('ASCL1 x hiNPCs_4w','ASCL1 x hiNPCs_6w','ASCL1 x NEUROG2'),]
EXC_DEGs_unique <- EXC_DEGs[!duplicated(EXC_DEGs$genes),]

###ASCL1 targets regulon network

regulons_targets <- read.table('results/regulons/regulons_EXC.txt') ## target genes of significant regulons

EXC_regulons <- read.csv("results/regulon_EXC_dataset.csv",header = T, stringsAsFactors = F,row.names = 1)

adj_EXC_sub <- adj_EXC[adj_EXC$TF %in% c('ASCL1',colnames(EXC_regulons)),] ### ASCL1 and regulons subset
adj_EXC_sub <- adj_EXC_sub[adj_EXC_sub$target %in% c(regulons_targets$V1,colnames(EXC_regulons)),]
adj_EXC_sub <- adj_EXC_sub[adj_EXC_sub$target %in% EXC_DEGs_unique$genes,]
adj_EXC_sub <- adj_EXC_sub[adj_EXC_sub$importance > 1,]

colnames(adj_EXC_sub) <- c('regulon','genes','enrichment')


motifs <- unique(as.character(adj_EXC_sub[,1]))
genes <- unique(as.character(adj_EXC_sub[,2]))

nodes <- data.frame(genes=c(motifs, genes),   
                    type=c(rep("Regulon", length(motifs)), rep("Target", length(genes))))
nodes <- nodes[!duplicated(nodes$genes),]


EXC_DEGs_unique <- EXC_DEGs_unique[EXC_DEGs_unique$genes %in% nodes$genes,]
nodes <- merge(nodes,EXC_DEGs_unique,by='genes',)
nodes <- nodes[!duplicated(nodes$genes),]

### Network creation
tal <- igraph::graph_from_data_frame(d = adj_EXC_sub,directed = T)

tal <- att.mapv(g = tal,dat=nodes,refcol = 1)

tal <- att.adde(tal, "edgeWidth", value = 1)

tal <- att.setv(g = tal,from = 'comparasion',to='nodeColor')

nodes_processed <- igraph::as_data_frame(tal, what = 'vertices')
nodes_processed$type <- ifelse(is.na(nodes_processed$type),'Regulon','Target')
nodes_processed$type <- ifelse(nodes_processed$name %in% c('ASCL1',colnames(EXC_regulons)),'Regulon','Target')
nodes_processed <- nodes_processed %>% mutate(color = case_when(
  nodeColor == '#B3B3B3' ~ '#5E39A3',
  TRUE ~ nodeColor
))


### Network customization 

tal <- igraph::graph_from_data_frame(d = adj_EXC_sub,directed = T)

tal <- att.mapv(g = tal,dat=nodes_processed,refcol = 1)

tal <- att.adde(tal, "edgeWidth", value = 2)
tal <- att.addv(tal,'nodeSize',value = 70)

tal <- att.setv(g = tal,from = 'comparasion',to='nodeColor',cols = c('#00008B','#7FD27F','#8B0000'))

tal <- att.setv(g = tal,from = 'type',to='nodeShape',shapes = c('DIAMOND','ELLIPSE'))
rdp <- RedPort()
calld(rdp)
resetd(rdp)

addGraph(rdp, tal, gcoord=c(10,25), gscale=20, theme='tm1', gzoom=30)

scl <- tal$legNodeColor$scale
leg <- tal$legNodeColor$legend 
addLegend.color(rdp, colvec=scl, labvec=leg, title="Comparasion")


####

logo_motiff_EXC <- readr::read_csv("~/GRN/results/logo_motiff_EXC.tsv")




