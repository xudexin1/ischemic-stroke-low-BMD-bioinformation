rm(list = ls())
setwd('/home/data/project/65.CCPJ-2003-6/')
file.dir = '07_roc'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
#IS 训练集------
library(magrittr)
library(ggplot2)
dat <- read.csv('../00_rawdata/dat.GSE58294.csv',row.names = 1)
group <- read.csv('../00_rawdata/group.GSE58294.csv')
table(group$group)
dat <- dat[,group$sample]
hubgene <- read.csv('../06_hubgene/hubgene.csv')
group$group = factor(group$group, levels = c("IS", "Control"))
hub_exp<-dat[hubgene$symbol,]
library(ggplot2)
library(pROC)
library(ROCR)
roc_list <- list()
for (gene in rownames(hub_exp)) {
  exp_vec <- as.numeric(hub_exp[gene, ])
  response <- as.numeric(group$group == "IS")  
  
  roc_obj <- pROC::roc(response, exp_vec, quiet = TRUE)
  roc_list[[gene]] <- roc_obj
}
library(IOBR)
colors <- c("#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pdf(file = '01.IS_GSE58294_roc.pdf',width = 5,height = 5)
plot(roc_list[[1]], col = colors[1], lwd = 2, main = "IS Corhort: ROC curves - GSE58294", legacy.axes = TRUE)

for (i in 2:length(roc_list)) {
  plot(roc_list[[i]], col = colors[i], lwd = 2, add = TRUE)
}

legend_text <- sapply(names(roc_list), function(gene) {
  auc_val <- round(auc(roc_list[[gene]]), 3)
  paste0(gene, " (AUC=", auc_val, ")")
})
legend("bottomright", legend = legend_text, col = colors, lwd = 2, cex = 0.8)
dev.off()

png(file = '01.IS_GSE58294_roc.png',width = 5,height = 5,units = 'in',res = 300,bg = 'white')
plot(roc_list[[1]], col = colors[1], lwd = 2, main = "IS Corhort: ROC curves - GSE58294", legacy.axes = TRUE)

for (i in 2:length(roc_list)) {
  plot(roc_list[[i]], col = colors[i], lwd = 2, add = TRUE)
}

legend_text <- sapply(names(roc_list), function(gene) {
  auc_val <- round(auc(roc_list[[gene]]), 3)
  paste0(gene, " (AUC=", auc_val, ")")
})
legend("bottomright", legend = legend_text, col = colors, lwd = 2, cex = 0.8)
dev.off()
# 
#OP 训练集------
library(magrittr)
library(ggplot2)
dat <- read.csv('../00_rawdata/dat.GSE56814.csv',row.names = 1)
group <- read.csv('../00_rawdata/group.GSE56814.csv')
table(group$group)
dat <- dat[,group$sample]
hubgene <- read.csv('../06_hubgene/hubgene.csv')

group$group = factor(group$group, levels = c("low BMD", "Control"),)
hub_exp<-dat[hubgene$symbol,]
library(ggplot2)
library(pROC)
library(ROCR)
roc_list <- list()
for (gene in rownames(hub_exp)) {
  exp_vec <- as.numeric(hub_exp[gene, ])
  response <- as.numeric(group$group == "low BMD")  
  
  roc_obj <- pROC::roc(response, exp_vec, quiet = TRUE)
  roc_list[[gene]] <- roc_obj
}
library(IOBR)
colors <- c("#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pdf(file = '02.OP_GSE56814_roc.pdf',width = 5,height = 5)
plot(roc_list[[1]], col = colors[1], lwd = 2, main = "low BMD Corhort: ROC curves - GSE56814", legacy.axes = TRUE)
for (i in 2:length(roc_list)) {
  plot(roc_list[[i]], col = colors[i], lwd = 2, add = TRUE)
}
legend_text <- sapply(names(roc_list), function(gene) {
  auc_val <- round(auc(roc_list[[gene]]), 3)
  paste0(gene, " (AUC=", auc_val, ")")
})
legend("bottomright", legend = legend_text, col = colors, lwd = 2, cex = 0.8)
dev.off()

png(file = '02.OP_GSE56814_roc.png',width = 5,height = 5,units = 'in',res = 300,bg = 'white')
plot(roc_list[[1]], col = colors[1], lwd = 2, main = "low BMD Corhort: ROC curves - GSE56814", legacy.axes = TRUE)

for (i in 2:length(roc_list)) {
  plot(roc_list[[i]], col = colors[i], lwd = 2, add = TRUE)
}

legend_text <- sapply(names(roc_list), function(gene) {
  auc_val <- round(auc(roc_list[[gene]]), 3)
  paste0(gene, " (AUC=", auc_val, ")")
})
legend("bottomright", legend = legend_text, col = colors, lwd = 2, cex = 0.8)
dev.off()
# 
#IS 验证集------
library(magrittr)
library(ggplot2)
dat <- read.csv('../00_rawdata/dat.GSE16561.csv',row.names = 1)
dat['CENPN',] <- dat['C16ORF60',]
group <- read.csv('../00_rawdata/group.GSE16561.csv')
table(group$group)
dat <- dat[,group$sample]
hubgene <- read.csv('../06_hubgene/hubgene.csv')
group$group = factor(group$group, levels = c("IS", "Control"))
hub_exp<-dat[hubgene$symbol,]
hub_exp<-hub_exp[-6,]
library(ggplot2)
library(pROC)
library(ROCR)
roc_list <- list()

for (gene in rownames(hub_exp)) {
  
  exp_vec <- as.numeric(hub_exp[gene, ])
  response <- as.numeric(group$group == "IS")
  roc_obj <- pROC::roc(
    response,
    exp_vec,
    quiet = TRUE,
    direction = "<"
  )
  
  if (as.numeric(pROC::auc(roc_obj)) < 0.5) {
    roc_obj <- pROC::roc(
      response,
      exp_vec,
      quiet = TRUE,
      direction = ">"
    )
  }
  
  roc_list[[gene]] <- roc_obj
}
library(IOBR)

colors <- c(
  "#E69F00", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7"
)

pdf(file = "03.IS_GSE16561_roc.pdf", width = 5, height = 5)
plot(
  roc_list[[1]],
  col = colors[1],
  lwd = 2,
  main = "IS Cohort: ROC curves - GSE16561",
  legacy.axes = TRUE
)

for (i in 2:length(roc_list)) {
  plot(
    roc_list[[i]],
    col = colors[i],
    lwd = 2,
    add = TRUE
  )
}

legend_text <- sapply(names(roc_list), function(gene) {
  auc_val <- round(as.numeric(pROC::auc(roc_list[[gene]])), 3)
  paste0(gene, " (AUC=", auc_val, ")")
})

legend(
  "bottomright",
  legend = legend_text,
  col = colors[seq_along(roc_list)],
  lwd = 2,
  cex = 0.8
)

dev.off()

png(file = '03.IS_GSE16561_roc.png',width = 5,height = 5,units = 'in',res = 300,bg = 'white')
plot(
  roc_list[[1]],
  col = colors[1],
  lwd = 2,
  main = "IS Cohort: ROC curves - GSE16561",
  legacy.axes = TRUE
)

for (i in 2:length(roc_list)) {
  plot(
    roc_list[[i]],
    col = colors[i],
    lwd = 2,
    add = TRUE
  )
}

legend_text <- sapply(names(roc_list), function(gene) {
  auc_val <- round(as.numeric(pROC::auc(roc_list[[gene]])), 3)
  paste0(gene, " (AUC=", auc_val, ")")
})

legend(
  "bottomright",
  legend = legend_text,
  col = colors[seq_along(roc_list)],
  lwd = 2,
  cex = 0.8
)

dev.off()
# 
#OP 验证集------
library(magrittr)
library(ggplot2)
dat <- read.csv('../00_rawdata/dat.GSE56815.csv',row.names = 1)
group <- read.csv('../00_rawdata/group.GSE56815.csv')
table(group$group)
dat <- dat[,group$sample]
hubgene <- read.csv('../06_hubgene/hubgene.csv')
group$group = factor(group$group, levels = c("low BMD", "Control"))
hub_exp<-dat[hubgene$symbol,]
library(ggplot2)
library(pROC)
library(ROCR)
roc_list <- list()
for (gene in rownames(hub_exp)) {
  exp_vec <- as.numeric(hub_exp[gene, ])
  response <- as.numeric(group$group == "low BMD")  
  
  roc_obj <- pROC::roc(response, exp_vec, quiet = TRUE)
  roc_list[[gene]] <- roc_obj
}
library(IOBR)
colors <- c("#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pdf(file = '04.OP_GSE56815_roc.pdf',width = 5,height = 5)
plot(roc_list[[1]], col = colors[1], lwd = 2, main = "low BMD Corhort: ROC curves - GSE56815", legacy.axes = TRUE)

for (i in 2:length(roc_list)) {
  plot(roc_list[[i]], col = colors[i], lwd = 2, add = TRUE)
}

legend_text <- sapply(names(roc_list), function(gene) {
  auc_val <- round(auc(roc_list[[gene]]), 3)
  paste0(gene, " (AUC=", auc_val, ")")
})
legend("bottomright", legend = legend_text, col = colors, lwd = 2, cex = 0.8)
dev.off()

png(file = '04.OP_GSE56815_roc.png',width = 5,height = 5,units = 'in',res = 300,bg = 'white')
plot(roc_list[[1]], col = colors[1], lwd = 2, main = "low BMD Corhort: ROC curves - GSE56815", legacy.axes = TRUE)
for (i in 2:length(roc_list)) {
  plot(roc_list[[i]], col = colors[i], lwd = 2, add = TRUE)
}

legend_text <- sapply(names(roc_list), function(gene) {
  auc_val <- round(auc(roc_list[[gene]]), 3)
  paste0(gene, " (AUC=", auc_val, ")")
})
legend("bottomright", legend = legend_text, col = colors, lwd = 2, cex = 0.8)
dev.off()
# 
