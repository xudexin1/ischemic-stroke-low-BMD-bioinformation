rm(list = ls())
setwd('/home/data/project/65.CCPJ-2003-6/')
file.dir = '03_enrichment'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
##IS-------
deg <- read.csv('../01_DEGs_IS/DEG_sig.csv',row.names = 1)
gene_transform <- bitr(rownames(deg),
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 1,
                readable = TRUE)
write.csv(ego,file = "IS/GO.csv",row.names = F)
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 70) +
  facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')+scale_fill_gradient(low = '#BB5566',high = '#80B1D3')
go_bar
ggsave('IS/01.GO_bar.png',go_bar,width =9,height = 6,units = 'in',dpi = 300,bg = 'white')
ggsave('IS/01.GO_bar.pdf',go_bar,width =9,height = 6)
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(kk@result,file = "IS/KEGG.csv",row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot

library(ggrepel)
kk_result <- kk@result[which(kk@result$p.adjust < 0.05),]
p <- ggplot(kk_result,aes(x=-log10(p.adjust),y=Description))+
  geom_point(aes(size=Count,color=Description),alpha=2)+
  scale_size(range = c(10,25))+
  theme_bw()+
  theme(legend.position = c('none'))+
  geom_text_repel(aes(label=Description),size=3,segment.color='black',show.legend = T)+
  theme(axis.text.y = element_blank())  +
  scale_color_manual(values = c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
    "#0072B2", "#D55E00", "#CC79A7", "#999999"
  )) 
p   
ggsave(filename = 'IS/02.KEGG_dot.pdf',p,w=6,h=5)
ggsave(filename = 'IS/02.KEGG_dot.png',p,w=6,h=5,dpi=300,units = 'in')


##OP-------
deg <- read.csv('../02_DEGs_OP/DEG_sig.csv',row.names = 1)
gene_transform <- bitr(rownames(deg),
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "none",
                pvalueCutoff = 0.05,
                qvalueCutoff = 1,
                readable = TRUE)
write.csv(ego,file = "OP/GO.csv",row.names = F)
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 70,color = 'pvalue') +
  facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')+scale_fill_gradient(low = '#BB5566',high = '#80B1D3')
go_bar
ggsave('OP/01.GO_bar.png',go_bar,width =9,height = 6,units = 'in',dpi = 300,bg = 'white')
ggsave('OP/01.GO_bar.pdf',go_bar,width =9,height = 6)
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "none",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(kk@result,file = "OP/KEGG.csv",row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot

library(ggrepel)
kk_result <- kk@result[which(kk@result$p.adjust < 0.05),]
kk_result <- kk_result[c(1:10),]
p <- ggplot(kk_result,aes(x=-log10(p.adjust),y=Description))+
  geom_point(aes(size=Count,color=Description),alpha=2)+
  scale_size(range = c(10,25))+
  theme_bw()+
  theme(legend.position = c('none'))+
  geom_text_repel(aes(label=Description),size=3,segment.color='black',show.legend = T)+
  theme(axis.text.y = element_blank())  +
  scale_color_manual(values = c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
    "#0072B2", "#D55E00", "#CC79A7", "#999999",'#332288',"#661100"
  )) 
p   
ggsave(filename = 'OP/02.KEGG_dot.pdf',p,w=6,h=5)
ggsave(filename = 'OP/02.KEGG_dot.png',p,w=6,h=5,dpi=300,units = 'in')
