rm(list = ls())
setwd('/home/data/project/65.CCPJ-2003-6/')
file.dir = '13_expression'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(tidyverse)
## 训练集 IS-----
dat <- read.csv('../00_rawdata/dat.GSE58294.csv',row.names = 1,check.names = F)
group <- read.csv('../00_rawdata/group.GSE58294.csv')
table(group$group)
hubgene <- read.csv('../06_hubgene/hubgene.csv')
control.sample <- group$sample[which(group$group=='Control')]
hub_exp<-dat[hubgene$symbol,]
hub_exp2 <- as.data.frame(hub_exp) #%>% na.omit()
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'Control','IS')

hub_exp2$Group <- factor(hub_exp2$Group,levels = c('Control','IS'))
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
stat.test$p<-ifelse(stat.test$p<0.0001,"****",
                    ifelse(stat.test$p<0.001,"***",
                           ifelse(stat.test$p<0.01,"**",
                                  ifelse(stat.test$p<0.05,"*",'ns'))))
colnames(hub_exp2)
exp_plot <- ggplot(hub_exp2, aes(x = Group, 
                                 y = expr, 
                                 fill = Group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15,
               position = position_dodge(0.75)) +
  
  geom_boxplot(width = 0.5,
               position = position_dodge(0.75),
               outlier.shape = NA,
               color = "black", 
               lwd = 0.3) +
  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
              size = 0.8,
              shape = 21,
              stroke = 0.2,
              color = "black") +
  
  scale_fill_manual(values = c("#56B4E9", "#009E73")) +  
  stat_compare_means(aes(group = Group), 
                     label = "p.signif", method = 'wilcox.test',
                     hide.ns = F, 
                     size = 4, 
                     label.y.npc = "top",label.x = 1.4) +
  
  labs(title = "IS Cohort: Expression Level in GSE58294", x = "", y = "Expression", size = 20) + 
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5,size = 16),
    axis.text.x = element_text(angle = 75, hjust = 1,  size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey60", size = 0.4)
  )+facet_wrap(~Symbol,scales = 'free',nrow = 1)
exp_plot
ggsave('01.expression_GSE58294_IS.pdf',exp_plot,width = 8,height = 5)
ggsave('01.expression_GSE58294_IS.png',exp_plot,width = 8,height = 5,units = 'in',dpi = 300)

## 验证集IS-----
dat <- read.csv('../00_rawdata/dat.GSE16561.csv',row.names = 1,check.names = F)
range(dat)
dat['CENPN',] <- dat['C16ORF60',]
allgene <- data.frame(rownames(dat))
group <- read.csv('../00_rawdata/group.GSE16561.csv')
table(group$group)
group$sample
control.sample <- group$sample[which(group$group=='Control')]
hubgene <- read.csv('../06_hubgene/hubgene.csv')

hub_exp<-dat[hubgene$symbol,]

hub_exp2 <- as.data.frame(hub_exp) %>% na.omit()
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'Control','IS')

hub_exp2$Group <- factor(hub_exp2$Group,levels = c('Control','IS'))
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
stat.test$p<-ifelse(stat.test$p<0.0001,"****",
                    ifelse(stat.test$p<0.001,"***",
                           ifelse(stat.test$p<0.01,"**",
                                  ifelse(stat.test$p<0.05,"*",'ns'))))
colnames(hub_exp2)
exp_plot <- ggplot(hub_exp2, aes(x = Group, 
                                 y = expr, 
                                 fill = Group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15,
               position = position_dodge(0.75)) +
  
  geom_boxplot(width = 0.5,
               position = position_dodge(0.75),
               outlier.shape = NA,
               color = "black", 
               lwd = 0.3) +
  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
              size = 0.8,
              shape = 21,
              stroke = 0.2,
              color = "black") +
  
  scale_fill_manual(values = c("#F0E442", "#D55E00")) +  
  stat_compare_means(aes(group = Group), 
                     label = "p.signif", method = 'wilcox.test',
                     hide.ns = F, 
                     size = 4, 
                     label.y.npc = "top",label.x = 1.4) +
  
  # 标签设置
  labs(title = "IS Cohort: Expression Level in GSE16561", x = "", y = "Expression", size = 20) + 
  
  # 美化主题
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5,size = 16),
    axis.text.x = element_text(angle = 75, hjust = 1,  size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey60", size = 0.4)
  )+facet_wrap(~Symbol,scales = 'free',nrow = 1)
exp_plot
ggsave('02.expression_GSE16561_IS.pdf',exp_plot,width = 8,height = 5)
ggsave('02.expression_GSE16561_IS.png',exp_plot,width = 8,height = 5,units = 'in',dpi = 300)


##训练集OP------
dat <- read.csv('../00_rawdata/dat.GSE56814.csv',row.names = 1,check.names = F)
group <- read.csv('../00_rawdata/group.GSE56814.csv')
table(group$group)
control.sample <- group$sample[which(group$group=='Control')]
hub_exp<-dat[hubgene$symbol,]
hub_exp2 <- as.data.frame(hub_exp) %>% na.omit()
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'Control','OP')

hub_exp2$Group <- factor(hub_exp2$Group,levels = c('Control','OP'))
table(hub_exp2$Group)
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
stat.test$p<-ifelse(stat.test$p<0.0001,"****",
                    ifelse(stat.test$p<0.001,"***",
                           ifelse(stat.test$p<0.01,"**",
                                  ifelse(stat.test$p<0.05,"*",'ns'))))
colnames(hub_exp2)
exp_plot <- ggplot(hub_exp2, aes(x = Group, 
                                 y = expr, 
                                 fill = Group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15,
               position = position_dodge(0.75)) +
  
  geom_boxplot(width = 0.5,
               position = position_dodge(0.75),
               outlier.shape = NA,
               color = "black", 
               lwd = 0.3) +
  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
              size = 0.8,
              shape = 21,
              stroke = 0.2,
              color = "black") +
  
  scale_fill_manual(values = c("#56B4E9", "#009E73")) + 
  stat_compare_means(aes(group = Group), 
                     label = "p.signif", method = 't.test',
                     hide.ns = F, 
                     size = 4, 
                     label.y.npc = "top",label.x = 1.4) +
  
  labs(title = "OP Cohort: Expression Level in GSE56814", x = "", y = "Expression", size = 20) + 
  
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5,size = 16),
    axis.text.x = element_text(angle = 75, hjust = 1,  size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey60", size = 0.4)
  )+facet_wrap(~Symbol,scales = 'free',nrow = 1)
exp_plot
ggsave('03.expression_GSE56814_OP.pdf',exp_plot,width = 8,height = 5)
ggsave('03.expression_GSE56814_OP.png',exp_plot,width = 8,height = 5,units = 'in',dpi = 300)

## OP test-----
dat <- read.csv('../00_rawdata/dat.GSE56815.csv',row.names = 1,check.names = F)
allgene <- data.frame(rownames(dat))
group <- read.csv('../00_rawdata/group.GSE56815.csv')
table(group$group)
control.sample <- group$sample[which(group$group=='Control')]
hub_exp<-dat[hubgene$symbol,]
hub_exp2 <- as.data.frame(hub_exp) %>% na.omit()
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'Control','OP')

hub_exp2$Group <- factor(hub_exp2$Group,levels = c('Control','OP'))
table(hub_exp2$Group)
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
stat.test$p<-ifelse(stat.test$p<0.0001,"****",
                    ifelse(stat.test$p<0.001,"***",
                           ifelse(stat.test$p<0.01,"**",
                                  ifelse(stat.test$p<0.05,"*",'ns'))))
colnames(hub_exp2)
exp_plot <- ggplot(hub_exp2, aes(x = Group, 
                                 y = expr, 
                                 fill = Group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.15,
               position = position_dodge(0.75)) +
  
  geom_boxplot(width = 0.5,
               position = position_dodge(0.75),
               outlier.shape = NA,
               color = "black", 
               lwd = 0.3) +
  
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
              size = 0.8,
              shape = 21,
              stroke = 0.2,
              color = "black") +
  
  scale_fill_manual(values = c("#F0E442", "#D55E00")) +
  stat_compare_means(aes(group = Group), 
                     label = "p.signif", method = 'wilcox.test',
                     hide.ns = F, 
                     size = 4, 
                     label.y.npc = "top",label.x = 1.4) +
  
  labs(title = "Expression Level in GSE56815", x = "", y = "Expression", size = 20) + 
  
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5,size = 16),
    axis.text.x = element_text(angle = 75, hjust = 1,  size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey60", size = 0.4)
  )+facet_wrap(~Symbol,scales = 'free',nrow = 1)
exp_plot
ggsave('04.expression_GSE56815_OP.pdf',exp_plot,width = 8,height = 5)
ggsave('04.expression_GSE56815_OP.png',exp_plot,width = 8,height = 5,units = 'in',dpi = 300)
