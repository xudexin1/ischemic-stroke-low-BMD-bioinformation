rm(list = ls())
setwd('/home/data/project/65.CCPJ-2003-6/')
file.dir = '00_rawdata'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
# IS GSE58294--------
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO("GSE58294",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL570",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(a)
group <- data.frame(sample = pd$geo_accession,group = pd$`group:ch1`)
table(group$group)
group$group <- ifelse(group$group == 'Cardioembolic Stroke','IS','Control')
group <- group[order(group$group),]
dat <- dat[,group$sample]
range(dat)
write.csv(dat,file = 'dat.GSE58294.csv',row.names = T)
write.csv(group,file = 'group.GSE58294.csv',row.names = F)

#IS-2 GSE16561------
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO("GSE16561",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr <- read.delim2('GSE16561_RAW.txt.gz',row.names = 1)
gpl<-getGEO("GPL6883",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','ILMN_Gene')%>%
  filter('ILMN_Gene'!='')%>%
  separate('ILMN_Gene',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(a)
group <- data.frame(sample = colnames(dat),group = colnames(dat))
table(group$group)
group$group <- ifelse(grepl('Stroke',group$group),'IS','Control')
# group$group <- ifelse(group$group == 'Stroke','IS','Control')
group <- group[order(group$group),]
dat <- dat[,group$sample]
write.csv(dat,file = 'dat.GSE16561.csv',row.names = T)
write.csv(group,file = 'group.GSE16561.csv',row.names = F)


#IS-3 GSE9877 ------
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO("GSE9877",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL96",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(a)
group <- data.frame(sample = pd$geo_accession,group = pd$title)
table(group$group)
group$group <- ifelse(grepl('normal',group$group),'Control','IS')
group <- group[order(group$group),]
dat <- dat[,group$sample]
range(dat)
dat <- log2(dat + 1)
write.csv(dat,file = 'dat.GSE9877.csv',row.names = T)
write.csv(group,file = 'group.GSE9877.csv',row.names = F)

#OP-1 GSE56815------
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO("GSE56815",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL96",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(a)
group <- data.frame(sample = pd$geo_accession,group = pd$`bone mineral density:ch1`)
table(group$group)
group$group <- ifelse(group$group == 'low BMD','low BMD','Control')
group <- group[order(group$group),]
dat <- dat[,group$sample]
range(dat)
write.csv(dat,file = 'dat.GSE56815.csv',row.names = T)
write.csv(group,file = 'group.GSE56815.csv',row.names = F)


#OP-2 GSE56814-----
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO("GSE56814",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL5175",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','gene_assignment')%>%
  filter('gene_assignment'!='')%>%
  separate('gene_assignment',c('drop','symbol'),sep = ' // ')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',] %>% na.omit()
allid <- data.frame(rownames(expr))
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(a)
group <- data.frame(sample = pd$geo_accession,group = pd$`bone mineral density:ch1`)
table(group$group)
group$group <- ifelse(group$group == 'OP','low BMD','Control')
group <- group[order(group$group),]
dat <- dat[,group$sample]
range(dat)
write.csv(dat,file = 'dat.GSE56814.csv',row.names = T)
write.csv(group,file = 'group.GSE56814.csv',row.names = F)

