rm(list = ls())
setwd('/home/data/project/65.CCPJ-2003-6/')
file.dir = '12_drug'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
drug <- read_tsv('gene_interaction_results.tsv')
table(drug$gene)
drug <- drug[which(drug$`interaction score` > 0.5),]
write.csv(drug,file = 'drug.select.csv',row.names = F)
length(unique(drug$drug))
colnames(drug)
str(drug)
edges <- drug[, c("gene", "drug")]
colnames(edges) <- c("from", "to")
graph <- graph_from_data_frame(edges, directed = TRUE)
deg <- degree(graph, mode = 'all')  # 获取度数
set.seed(123)
p <- ggraph(graph, layout = 'fr') +
  geom_edge_link(arrow = arrow(length = unit(4, 'mm')),
                 end_cap = circle(3, 'mm'),
                 edge_alpha = 0.5) +
  geom_node_point(aes(color = deg), size = 4) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3,max.overlaps = Inf) +
  scale_color_viridis_c(name = "Degree", option = "D", direction = -1) +
  theme_void() +
  labs(title = "mRNA-Drug Network")
p
ggsave(filename = '01.mRNA-drug.pdf',p,width = 6,height = 5)
ggsave(filename = '01.mRNA-drug.png',p,width = 6,height = 5,units = 'in',dpi = 300,bg = 'white')



library(igraph)
library(ggraph)
library(tidyverse)
edges <- drug[, c("gene", "drug")]
colnames(edges) <- c("from", "to")

genes <- unique(drug$gene)
drugs <- unique(drug$drug)

nodes <- data.frame(
  name = unique(c(genes, drugs)),
  type = ifelse(unique(c(genes, drugs)) %in% genes, "Gene", "Drug")
)

drug_info <- drug %>%
  select(drug, `regulatory approval`) %>%
  distinct() %>%
  rename(name = drug)

nodes <- left_join(nodes, drug_info, by = "name")

nodes$approval_status <- ifelse(nodes$type == "Drug",
                                nodes$`regulatory approval`,
                                "Gene")

graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
deg <- degree(graph, mode = "all")
V(graph)$deg <- deg
V(graph)$type <- V(graph)$type
V(graph)$approval_status <- V(graph)$approval_status

set.seed(123)
p <- ggraph(graph, layout = 'fr') +
  geom_edge_link(arrow = arrow(length = unit(3, 'mm')),
                 end_cap = circle(3, 'mm'),
                 edge_alpha = 0.3) +
  geom_node_point(aes(color = approval_status, shape = type, size = deg)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = Inf) +
  scale_color_manual(
    values = c(
      "Approved" = "#1b9e77",
      "Not Approved" = "#d95f02",
      "Gene" = "#7570b3"
    ),
    breaks = c("Approved", "Not Approved"),  # 仅这两类显示在图例
    name = "Approval"
  ) +
  scale_shape_manual(values = c("Gene" = 16, "Drug" = 15)) +
  scale_size(range = c(3, 8)) +
  theme_void() +
  labs(title = "mRNA-Drug Network",
       shape = "Node Type",
       size = "Degree")

p
ggsave(filename = "01.mRNA-drug-network.pdf", p, width = 8, height = 6)
ggsave(filename = "01.mRNA-drug-network.png", p, width = 8, height = 6, dpi = 300, bg = "white")
