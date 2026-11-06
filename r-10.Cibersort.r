rm(list = ls())
setwd('/home/data/project/65.CCPJ-2003-6/')
file.dir = '10_Cibersort'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

library(tidyverse)
dat <- read.csv('../00_rawdata/dat.GSE58294.csv',row.names = 1,check.names = F)
group <- read.csv('../00_rawdata/group.GSE58294.csv')
str(group)
dat <- dat[,group$sample]
library(tidyverse)
library(immunedeconv)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
set_cibersort_binary("../../../files/TIIC/CIBERSORT.R")
set_cibersort_mat("../../../files/TIIC/LM22.txt")
res.cibersort <- deconvolute(as.matrix(dat), method = 'cibersort')
save(file = 'IS/cibersort.Rdata', res.cibersort)

mypalette <- colorRampPalette(brewer.pal(9, "Set1"))
display.brewer.all()
plot_cibersort <- function(data, group) {
  data %>%
    gather(sample, fraction, -cell_type) %>%
    merge(group, by = 'sample') %>%
    ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
    geom_bar(position = 'stack', stat = 'identity') +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    labs(x = '', y = 'Relative Percent', fill = '') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'right') +
    scale_fill_manual(values = mypalette(22)) +
    facet_grid(~group, scales = "free", space = "free")+
    guides(fill = guide_legend(ncol = 1))
}

pdf('IS/01.cibersort.box.pdf', width = 9, height = 6)
print(plot_cibersort(res.cibersort, group))
dev.off()

png('IS/01.cibersort.box.png', width = 9, height = 6, units = 'in', res = 300,bg = 'white')
print(plot_cibersort(res.cibersort, group))
dev.off()
res.cibersort2 <- res.cibersort %>% column_to_rownames(var = 'cell_type')
keep <- rowSums(res.cibersort2 > 0) >= floor(0.10 * ncol(res.cibersort2))
res.cibersort2 <- res.cibersort2[keep, ]
dat.cibersort <- res.cibersort2 %>%
  t %>% as.data.frame() %>%
  rownames_to_column(var = "sample")
dat.cibersort <- merge(group, dat.cibersort, by = "sample")
dat.cibersort2 <- gather(dat.cibersort, ImmuneCell, Score, -c("sample", "group"))

stat_cibersort <- dat.cibersort2 %>%
  group_by(ImmuneCell) %>%
  wilcox_test(Score ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p")

write.csv(stat_cibersort, 'IS/stat.cibersort.csv', row.names = TRUE)
DE.cibersort <- stat_cibersort %>% filter(p < 0.05)
violin.cibersort <- dat.cibersort2 %>% filter(ImmuneCell %in% DE.cibersort$ImmuneCell)
violin.cibersort$group <- factor(violin.cibersort$group, levels = c('Control', 'IS'))
num_cells <- length(unique(violin.cibersort$ImmuneCell))
nrow_val <- ceiling(num_cells / 5)
cibersort_plot <- ggplot(violin.cibersort, aes(x = group, y = Score, fill = group)) +
  geom_violin(trim = FALSE, color = "black") +
  stat_boxplot(geom = "errorbar", width = 0.1, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA, fill = 'white') +
  scale_fill_manual(values = c('#E1C62F',"#A24057"), name = "Group") +
  labs(title = "Immune Cell", x = "", y = "Score") +
  stat_compare_means(label = "p.signif", method = 'wilcox.test', paired = FALSE, label.x = 1.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~ImmuneCell, scales = "free", nrow = nrow_val)
cibersort_plot
ggsave('IS/02.cibersort.plot.pdf', plot = cibersort_plot, width = 9, height = nrow_val*2+2)
ggsave('IS/02.cibersort.plot.png', plot = cibersort_plot, width = 10.5, height = nrow_val*2+3, dpi = 300)

##关键基因和关键免疫细胞相关性---------
library(ggcorrplot)
library(corrplot)
library(Hmisc)
hubgene <- read.csv('../06_hubgene/hubgene.csv')
hub.exp <- dat[hubgene$symbol,]
detiics <- res.cibersort2[which(rownames(res.cibersort2) %in% stat_cibersort$ImmuneCell[which(stat_cibersort$p < 0.05)]),]
nc<-t(rbind(hub.exp,detiics))
m=rcorr(nc,type = 'spearman')$r[1:nrow(hub.exp),(ncol(nc)-length(rownames(detiics))+1):ncol(nc)]
m<-t(m)

p=rcorr(nc,type = 'spearman')$P[1:nrow(hub.exp),(ncol(nc)-length(rownames(detiics))+1):ncol(nc)]
p<-t(p)
library(reshape2)
correlation <- reshape2::melt(m)
colnames(correlation) <- c('Cell','Gene','correlation')
pvalue <- reshape2::melt(p)
correlation$pvalue <- pvalue$value
write.csv(correlation,file = 'IS/correlation.csv',row.names = F)

library(reshape2)
library(ggplot2)

corr_long <- reshape2::melt(t(m))
pval_long <- reshape2::melt(t(p))
colnames(corr_long) <- c("Var1", "Var2", "cor")
colnames(pval_long) <- c("Var1", "Var2", "pval")
df <- merge(corr_long, pval_long, by = c("Var1", "Var2"))

df$sig <- cut(df$pval,
              breaks = c(-Inf,0.0001, 0.001, 0.01, 0.05, Inf),
              labels = c("****", "***", "**", "*", ""))

df$logp <- df$cor

pdf(file = 'IS/03.correlation.pdf',w=6,h=5)
ggplot(df, aes(x = Var1, y = Var2, fill = logp)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), size = 5, color = "black") +
  scale_fill_gradientn(colors = c("#377eb8", "white", "#e41a1c"),
                       name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank()) +
  labs(title = "IS Corhort: Spearman Correlation - GSE58294",
       x = NULL, y = NULL)#+coord_flip()

dev.off()


png(file = 'IS/03.correlation.png',w=6,h=5,units = 'in',res = 300,bg = 'white')
ggplot(df, aes(x = Var1, y = Var2, fill = logp)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), size = 5, color = "black") +
  scale_fill_gradientn(colors = c("#377eb8", "white", "#e41a1c"),
                       name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank()) +
  labs(title = "IS Corhort: Spearman Correlation - GSE58294",
       x = NULL, y = NULL)#+coord_flip()

dev.off()


##02 OP-cibersort------
dat <- read.csv('../00_rawdata/dat.GSE56814.csv',row.names = 1,check.names = F)
group <- read.csv('../00_rawdata/group.GSE56814.csv')
str(group)
dat <- dat[,group$sample]
library(tidyverse)
library(immunedeconv)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
set_cibersort_binary("../../../files/TIIC/CIBERSORT.R")
set_cibersort_mat("../../../files/TIIC/LM22.txt")
res.cibersort <- deconvolute(as.matrix(dat), method = 'cibersort')
save(file = 'OP/cibersort.Rdata', res.cibersort)

mypalette <- colorRampPalette(brewer.pal(9, "Set1"))
display.brewer.all()
plot_cibersort <- function(data, group) {
  data %>%
    gather(sample, fraction, -cell_type) %>%
    merge(group, by = 'sample') %>%
    ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
    geom_bar(position = 'stack', stat = 'identity') +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    labs(x = '', y = 'Relative Percent', fill = '') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'right') +
    scale_fill_manual(values = mypalette(22)) +
    facet_grid(~group, scales = "free", space = "free")+
    guides(fill = guide_legend(ncol = 1))
}

pdf('OP/01.cibersort.box.pdf', width = 9, height = 6)
print(plot_cibersort(res.cibersort, group))
dev.off()

png('OP/01.cibersort.box.png', width = 9, height = 6, units = 'in', res = 300,bg = 'white')
print(plot_cibersort(res.cibersort, group))
dev.off()

res.cibersort2 <- res.cibersort %>% column_to_rownames(var = 'cell_type')
keep <- rowSums(res.cibersort2 > 0) >= floor(0.10 * ncol(res.cibersort2))
res.cibersort2 <- res.cibersort2[keep, ]
dat.cibersort <- res.cibersort2 %>%
  t %>% as.data.frame() %>%
  rownames_to_column(var = "sample")
dat.cibersort <- merge(group, dat.cibersort, by = "sample")
dat.cibersort2 <- gather(dat.cibersort, ImmuneCell, Score, -c("sample", "group"))

stat_cibersort <- dat.cibersort2 %>%
  group_by(ImmuneCell) %>%
  wilcox_test(Score ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p")

write.csv(stat_cibersort, 'OP/stat.cibersort.csv', row.names = TRUE)
DE.cibersort <- stat_cibersort %>% filter(p < 0.05)
violin.cibersort <- dat.cibersort2 %>% filter(ImmuneCell %in% DE.cibersort$ImmuneCell)
violin.cibersort$group <- factor(violin.cibersort$group, levels = c('Control', 'OP'))
num_cells <- length(unique(violin.cibersort$ImmuneCell))
nrow_val <- ceiling(num_cells / 5)
cibersort_plot <- ggplot(violin.cibersort, aes(x = group, y = Score, fill = group)) +
  geom_violin(trim = FALSE, color = "black") +
  stat_boxplot(geom = "errorbar", width = 0.1, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA, fill = 'white') +
  scale_fill_manual(values = c('#E1C62F',"#A24057"), name = "Group") +
  labs(title = "Immune Cell", x = "", y = "Score") +
  stat_compare_means(label = "p.signif", method = 'wilcox.test', paired = FALSE, label.x = 1.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~ImmuneCell, scales = "free", nrow = nrow_val)
cibersort_plot
ggsave('OP/02.cibersort.plot.pdf', plot = cibersort_plot, width = 9, height = nrow_val*2+2)
ggsave('OP/02.cibersort.plot.png', plot = cibersort_plot, width = 10.5, height = nrow_val*2+3, dpi = 300)
library(ggcorrplot)
library(corrplot)
library(Hmisc)
hubgene <- read.csv('../06_hubgene/hubgene.csv')
hub.exp <- dat[hubgene$symbol,]
detiics <- res.cibersort2[which(rownames(res.cibersort2) %in% stat_cibersort$ImmuneCell[which(stat_cibersort$p < 0.05)]),]
nc<-t(rbind(hub.exp,detiics))
m=rcorr(nc,type = 'spearman')$r[1:nrow(hub.exp),(ncol(nc)-length(rownames(detiics))+1):ncol(nc)]
m<-t(m)

p=rcorr(nc,type = 'spearman')$P[1:nrow(hub.exp),(ncol(nc)-length(rownames(detiics))+1):ncol(nc)]
p<-t(p)
library(reshape2)
correlation <- reshape2::melt(m)
colnames(correlation) <- c('Cell','Gene','correlation')
pvalue <- reshape2::melt(p)
correlation$pvalue <- pvalue$value
write.csv(correlation,file = 'OP/correlation.csv',row.names = F)

library(reshape2)
library(ggplot2)

corr_long <- reshape2::melt(t(m))
pval_long <- reshape2::melt(t(p))
colnames(corr_long) <- c("Var1", "Var2", "cor")
colnames(pval_long) <- c("Var1", "Var2", "pval")

df <- merge(corr_long, pval_long, by = c("Var1", "Var2"))

df$sig <- cut(df$pval,
              breaks = c(-Inf,0.0001, 0.001, 0.01, 0.05, Inf),
              labels = c("****", "***", "**", "*", ""))
df$logp <- df$cor

pdf(file = 'OP/03.correlation.pdf',w=6,h=4)
ggplot(df, aes(x = Var1, y = Var2, fill = logp)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), size = 5, color = "black") +
  scale_fill_gradientn(colors = c("#377eb8", "white", "#e41a1c"),
                       name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank()) +
  labs(title = "OP Corhort: Spearman Correlation - GSE56814",
       x = NULL, y = NULL)#+coord_flip()

dev.off()


png(file = 'OP/03.correlation.png',w=6,h=4,units = 'in',res = 300,bg = 'white')
ggplot(df, aes(x = Var1, y = Var2, fill = logp)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), size = 5, color = "black") +
  scale_fill_gradientn(colors = c("#377eb8", "white", "#e41a1c"),
                       name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank()) +
  labs(title = "OP Corhort: Spearman Correlation - GSE56814",
       x = NULL, y = NULL)#+coord_flip()

dev.off()

