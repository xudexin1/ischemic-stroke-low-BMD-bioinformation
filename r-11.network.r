rm(list = ls())
setwd('/home/data/project/65.CCPJ-2003-6/')
file.dir = '11_network'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
gene2tf <- read.csv('gene2tf.csv')
length(unique(gene2tf$ID))
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
edges <- gene2tf[, c("ID", "Target")]
colnames(edges) <- c("from", "to")

# 构建网络
graph <- graph_from_data_frame(edges, directed = TRUE)
deg <- degree(graph, mode = 'all')  # 获取度数
set.seed(123)
# 绘图
p <- ggraph(graph, layout = 'fr') +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')),
                 end_cap = circle(3, 'mm'),
                 edge_alpha = 0.5) +
  geom_node_point(aes(color = deg), size = 4) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_viridis_c(name = "Degree", option = "G", direction = -1) +
  theme_void() +
  labs(title = "TF-mRNA Regulation Network")
p
ggsave(filename = '01.mRNA-TF.pdf',p,width = 6,height = 5)
ggsave(filename = '01.mRNA-TF.png',p,width = 6,height = 5,units = 'in',dpi = 300,bg = 'white')


##mirwalk-----
mirwalk <- read.csv('miRWalk_miRNA_Targets.csv')
mirwalk_validated <- mirwalk[mirwalk$validated != "", ]
colnames(mirwalk_validated)
edges <- mirwalk_validated[, c("mirnaid", "genesymbol")]
edges$link <- paste0(edges$mirnaid,edges$genesymbol)
edges <- edges[!duplicated(edges$link),] %>% select(-'link')
colnames(edges) <- c("from", "to")
length(unique(edges$from))
length(unique(edges$to))
graph <- graph_from_data_frame(edges, directed = TRUE)
deg <- degree(graph, mode = 'all')  # 获取度数
set.seed(123456)
p <- ggraph(graph, layout = 'fr') +
  geom_edge_link(arrow = arrow(length = unit(2, 'mm')),
                 end_cap = circle(3, 'mm'),
                 edge_alpha = 0.5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = Inf)+
  geom_node_point(aes(color = deg), size = 4) +
  scale_color_viridis_c(name = "Degree", option = "G", direction = -1) +
  theme_void() +
  labs(title = "miRAN-mRNA Regulation Network")
p
ggsave(filename = '02.miRNA-mRNA.pdf',p,width = 7,height = 6)

ggsave(filename = '02.miRNA-mRNA.png',p,width = 7,height = 6,units = 'in',dpi = 300,bg = 'white')
