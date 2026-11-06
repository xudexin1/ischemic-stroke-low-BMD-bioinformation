rm(list = ls())
setwd('/home/data/project/65.CCPJ-2003-6/')
file.dir = '08_GSEA'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(limma)
library(tidyverse)
library(magrittr)

##IS------
expr <- read.csv('../00_rawdata/dat.GSE58294.csv',row.names = 1)
gene <- read.csv('../06_hubgene/hubgene.csv')
gene <- gene$symbol
gene_expr <- expr[gene,] %>% t %>% as.data.frame
allgene_expr <- expr %>% t %>% as.data.frame

cor_res <- cor(gene_expr, allgene_expr,method = "spearman") %>% t %>% as.data.frame 
cor_res <- na.omit(cor_res)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

j <- 0
for (num in 1:length(gene)){
  i <- gene[num]
  print(i)
  j <- j+1
  dat <- data.frame(row.names = rownames(cor_res),
                    cor=cor_res[,i])
  data_sort <- dat %>% dplyr::arrange(desc(cor)) %>% na.omit()
  
  geneList <- data_sort$cor
  names(geneList) <- rownames(data_sort) 
  
  # kegg------------------------------------------------------------------
  kegg_set<- read.gmt("../../../files/GSVA/h.all.v2024.1.Hs.symbols.gmt")  
  set.seed(1)
  kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 1) 
  kegg_result <- kegg_gsea@result
  kegg_result <- kegg_result[kegg_result$p.adjust < 0.05,]
  kegg_result <- kegg_result[order(abs(kegg_result$enrichmentScore), decreasing = T),]
  write.csv(kegg_result,paste0('IS/','0',j,'.GSEA_Hallmark_result_',i,'.csv'),row.names = F)
  ifelse(length(rownames(kegg_result) )> 5,
         paths <- rownames(kegg_result[c(1:5),]),
         paths <- rownames(kegg_result)
  )
  my_colors <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")   
  my_colors <- my_colors[1:length(paths)]
  gsea <- gseaplot2(kegg_gsea,paths,color = my_colors,
                    title = i,
                    base_size = 11,
                    rel_heights = c(1.5, 0.3, 0.5))
  png(paste0('IS/','0',j,'.GSEA_Hallmark_',i,'.png'),w=8,h=6,units = "in",res = 300,bg = 'white')
  print(gsea)
  dev.off()
  pdf(paste0('IS/','0',j,'.GSEA_Hallmark_',i,'.pdf'),w=8,h=6)
  print(gsea)
  dev.off()
  
  num <- num+1
}


##OP------
expr <- read.csv('../00_rawdata/dat.GSE56814.csv',row.names = 1)
gene <- read.csv('../06_hubgene/hubgene.csv')
gene <- gene$symbol
gene_expr <- expr[gene,] %>% t %>% as.data.frame
allgene_expr <- expr %>% t %>% as.data.frame

cor_res <- cor(gene_expr, allgene_expr,method = "spearman") %>% t %>% as.data.frame 
cor_res <- na.omit(cor_res)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

j <- 0
for (num in 1:length(gene)){
  i <- gene[num]
  print(i)
  j <- j+1
  dat <- data.frame(row.names = rownames(cor_res),
                    cor=cor_res[,i])
  data_sort <- dat %>% dplyr::arrange(desc(cor)) %>% na.omit()
  
  geneList <- data_sort$cor
  names(geneList) <- rownames(data_sort) 
  
  # kegg------------------------------------------------------------------
  kegg_set<- read.gmt("../../../files/GSVA/h.all.v2024.1.Hs.symbols.gmt")  
  set.seed(1)
  kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 1)  
  kegg_result <- kegg_gsea@result
  kegg_result <- kegg_result[kegg_result$p.adjust < 0.05,]
  kegg_result <- kegg_result[order(abs(kegg_result$enrichmentScore), decreasing = T),]
  write.csv(kegg_result,paste0('OP/','0',j,'.GSEA_Hallmark_result_',i,'.csv'),row.names = F)
  ifelse(length(rownames(kegg_result) )> 5,
         paths <- rownames(kegg_result[c(1:5),]),
         paths <- rownames(kegg_result)
  )
  my_colors <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7")   
  my_colors <- my_colors[1:length(paths)]
  gsea <- gseaplot2(kegg_gsea,paths,color = my_colors,
                    title = i,
                    base_size = 11,
                    rel_heights = c(1.5, 0.3, 0.5))
  png(paste0('OP/','0',j,'.GSEA_Hallmark_',i,'.png'),w=8,h=6,units = "in",res = 300,bg = 'white')
  print(gsea)
  dev.off()
  pdf(paste0('OP/','0',j,'.GSEA_Hallmark_',i,'.pdf'),w=8,h=6)
  print(gsea)
  dev.off()
  
  num <- num+1
}

