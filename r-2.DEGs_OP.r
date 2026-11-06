rm(list = ls())
setwd('/home/data/project/65.CCPJ-2003-6/')
file.dir = '02_DEGs_OP'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(magrittr)
library(stringr)
library(limma)
df = read.csv('../00_rawdata/dat.GSE56814.csv',row.names = 1) 
range(df)
boxplot(df,las = 2)
df.group = read.csv('../00_rawdata/group.GSE56814.csv')
df = df[,df.group$sample]
table(df.group$group)
df.group$group <- gsub(' ','_',df.group$group)
df.group$group <- ifelse(df.group$group == 'Control','Control','low_BMD')

df.group$group = factor(df.group$group, levels = c("Control", "low_BMD"))
design.mat = cbind(Control = ifelse(df.group$group == "Control", 1, 0), 
                   low_BMD = ifelse(df.group$group == "Control", 0, 1))

contrast.mat = makeContrasts(contrasts="low_BMD-Control", levels=design.mat)
fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
DEG=na.omit(fit)
logFC_cutoff <- 0
DEG$change = as.factor(
  ifelse(DEG$P.Value <0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
summary(sig_diff$change) 
write.csv(DEG,file = "DEG_all.csv",row.names = T)
write.csv(sig_diff,file = "DEG_sig.csv",row.names = T)

library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)
library(ggrepel)
library(ggplot2)
library(ggrepel)
library(ggthemes)

dat_rep <- DEG[rownames(DEG) %in%
                 rownames(rbind(
                   head(sig_diff[order(sig_diff$logFC, decreasing = TRUE),], 10),
                   head(sig_diff[order(sig_diff$logFC, decreasing = FALSE),], 10)
                 )),]

volcano_plot <- ggplot(DEG, aes(x = logFC, y = -log10(P.Value), color = change)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("#80B1D3", "#e4dfd7", "#BB5566")) +
  # geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(data = dat_rep, aes(label = rownames(dat_rep)), 
                  size = 4, max.overlaps = 30,
                  box.padding = 0.5, point.padding = 0.4, 
                  segment.color = "gray30") +
  labs(x = "log2(Fold Change)", y = "-log10(adjusted P-value)",
       title = 'Differential Expression Volcano Plot:GSE56814') +
  theme_few(base_size = 14) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"))

volcano_plot

ggsave('01.volcano.png', volcano_plot, width = 8, height = 6, units = 'in', dpi = 300)
ggsave('01.volcano.pdf', volcano_plot, width = 8, height = 6)

