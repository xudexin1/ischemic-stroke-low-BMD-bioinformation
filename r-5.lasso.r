rm(list = ls())
setwd('/home/data/project/65.CCPJ-2003-6/')
file.dir = '05_Lasso'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
library(magrittr)
library(ggplot2)
library(caret)
library(glmnet)
# 训练集-----------
dat <- read.csv('../00_rawdata/dat.GSE58294.csv',row.names = 1)
dat <- dat %>% t() %>% as.data.frame()
candigene <- read.csv('../04_candigene/candigene.csv')
dat <- dat[,candigene$symbol]
group = read.csv('../00_rawdata/group.GSE58294.csv')
gene <- candigene
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
table(group$group)
dat$group<-factor(dat$group,levels = c('IS','Control'))

train_data <- dat
set.seed(123)
# 2. 训练Lasso回归模型
res.lasso <- cv.glmnet(as.matrix(train_data[-ncol(train_data)]), train_data$group, 
                       family = "binomial", type.measure = "deviance", nfolds = 10)

plot(res.lasso)
plot(res.lasso$glmnet.fit, xvar = 'lambda')
cv_data <- data.frame(
  log_lambda = log(res.lasso$lambda),
  cvm = res.lasso$cvm,
  cvup = res.lasso$cvup,
  cvlo = res.lasso$cvlo
)
cv_plot <- ggplot(cv_data, aes(x = log_lambda, y = cvm)) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), width = 0.2, color = "grey50") +
  geom_point(size = 2, color = "#104680") +
  geom_vline(xintercept = log(res.lasso$lambda.min), linetype = "dashed", color = "#ba2f7b") +
  geom_vline(xintercept = log(res.lasso$lambda.1se), linetype = "dashed", color = "#0f59a4") +
  labs(title = "IS Cohort: LASSO Cross Validation Curve - GSE58294",
       x = "log(Lambda)",
       y = "Mean Cross-Validated Error") +
  theme_minimal(base_size = 14)
cv_plot
ggsave("IS/01.lasso.CV.png", cv_plot, width = 8, height = 6, dpi = 300,bg = 'white')
ggsave("IS/01.lasso.CV.pdf", cv_plot, width = 8, height = 6)

library(broom)
df_coef <- as.data.frame(as.matrix(res.lasso$glmnet.fit$beta))
df_coef$gene <- rownames(df_coef)

library(reshape2)
df_long <- melt(df_coef, id.vars = "gene", variable.name = "Lambda_Index", value.name = "Coefficient")
df_long$Lambda <- rep(log(res.lasso$glmnet.fit$lambda), each = nrow(df_coef))
library(IOBR)
cols = palettes(category = "random",20,show_col=F)
coef_plot <- ggplot(df_long, aes(x = Lambda, y = Coefficient, color = gene)) +
  geom_line(size = 1) +
  geom_vline(xintercept = log(res.lasso$lambda.min), linetype = "dashed", color = "red") +
  labs(title = "IS Cohort: LASSO Coefficient Paths by Gene - GSE58294",
       x = "log(Lambda)",
       y = "Coefficient") +
  scale_color_manual(values = cols)+
  theme_minimal(base_size = 14) +
  theme(legend.position = "right") +
  guides(color = guide_legend(ncol = 2))
coef_plot

ggsave("IS/02.lasso.Coef.png", coef_plot, width = 8, height = 6, dpi = 300)
ggsave("IS/02.lasso.Coef.pdf", coef_plot, width = 8, height = 6)


l.coef<-coef(res.lasso$glmnet.fit,s=res.lasso$lambda.min,exact= F)
l.coef
coef.min = coef(res.lasso, s = "lambda.min")  
res.lasso$lambda.min
active.min = which(coef.min@i != 0)
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids <- lasso_geneids[-1]
lasso_geneids
write.csv(lasso_geneids,file = 'IS/lasso.gene.csv',row.names = F)
str(df_long)

beta_min <- as.matrix(coef(res.lasso, s = "lambda.min"))
beta_df <- data.frame(
  gene = rownames(beta_min),
  coef = as.numeric(beta_min)
)
beta_df <- beta_df[beta_df$coef != 0, ]
beta_df <- beta_df[-1,]
max_coef <- max(beta_df$coef) * 1.2

p <- ggplot(beta_df, aes(x = reorder(gene, coef), y = coef, fill = coef > 0)) +
  geom_col(width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#b5aa0f", "FALSE" = "#7ba39d")) +
  geom_label(aes(x = gene, y = max_coef, label = coef),
             hjust = 0, size = 3.5,
             fill = "#f7c173", color = "black", label.size = 0.3, label.padding = unit(0.15, "lines"))+
  
  labs(x = "Gene", y = "Coefficient", title = "IS Cohort: Gene Coefficient - GSE58294") +
  theme_bw() +
  theme(legend.position = "none") +
  
  scale_y_continuous(limits = c(min(beta_df$coef), max_coef * 1.5))
p

ggsave("IS/03.Coef.Bar.png", p, width = 8, height =5, dpi = 300,units = 'in',bg = 'white')
ggsave("IS/03.Coef.Bar.pdf", p, width = 8, height =5)

# OP-----------
dat <- read.csv('../00_rawdata/dat.GSE56814.csv',row.names = 1)
range(dat)
dat <- dat[candigene$symbol,] %>% t() %>% as.data.frame()
group = read.csv('../00_rawdata/group.GSE56814.csv')
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
table(group$group)
dat$group<-factor(dat$group,levels = c('low BMD','Control'))
train_data <- dat
set.seed(123)
# 2. 训练Lasso回归模型
res.lasso <- cv.glmnet(as.matrix(train_data[-ncol(train_data)]), train_data$group, 
                       family = "binomial", type.measure = "deviance", nfolds = 10)

plot(res.lasso)
plot(res.lasso$glmnet.fit, xvar = 'lambda')
cv_data <- data.frame(
  log_lambda = log(res.lasso$lambda),
  cvm = res.lasso$cvm,
  cvup = res.lasso$cvup,
  cvlo = res.lasso$cvlo
)
cv_plot <- ggplot(cv_data, aes(x = log_lambda, y = cvm)) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup), width = 0.2, color = "grey50") +
  geom_point(size = 2, color = "#104680") +
  geom_vline(xintercept = log(res.lasso$lambda.min), linetype = "dashed", color = "#ba2f7b") +
  geom_vline(xintercept = log(res.lasso$lambda.1se), linetype = "dashed", color = "#0f59a4") +
  labs(title = "low BMD Cohort: LASSO Cross Validation Curve - GSE56814",
       x = "log(Lambda)",
       y = "Mean Cross-Validated Error") +
  theme_minimal(base_size = 14)
cv_plot
ggsave("Low_BMD/01.lasso.CV.png", cv_plot, width = 8, height = 6, dpi = 300,bg = 'white')
ggsave("Low_BMD/01.lasso.CV.pdf", cv_plot, width = 8, height = 6)

library(broom)

df_coef <- as.data.frame(as.matrix(res.lasso$glmnet.fit$beta))
df_coef$gene <- rownames(df_coef)

library(reshape2)
df_long <- melt(df_coef, id.vars = "gene", variable.name = "Lambda_Index", value.name = "Coefficient")
df_long$Lambda <- rep(log(res.lasso$glmnet.fit$lambda), each = nrow(df_coef))
library(IOBR)
cols = palettes(category = "random",20,show_col=F)
coef_plot <- ggplot(df_long, aes(x = Lambda, y = Coefficient, color = gene)) +
  geom_line(size = 1) +
  geom_vline(xintercept = log(res.lasso$lambda.min), linetype = "dashed", color = "red") +
  labs(title = "low BMD Cohort: LASSO Coefficient Paths by Gene - GSE56814",
       x = "log(Lambda)",
       y = "Coefficient") +
  scale_color_manual(values = cols)+
  theme_minimal(base_size = 14) +
  theme(legend.position = "right") +
  guides(color = guide_legend(ncol = 2))
coef_plot

ggsave("Low_BMD/02.lasso.Coef.png", coef_plot, width = 8, height = 6, dpi = 300)
ggsave("Low_BMD/02.lasso.Coef.pdf", coef_plot, width = 8, height = 6)


l.coef<-coef(res.lasso$glmnet.fit,s=res.lasso$lambda.min,exact= F)
l.coef
coef.min = coef(res.lasso, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
res.lasso$lambda.min
active.min = which(coef.min@i != 0)
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids <- lasso_geneids[-1]
lasso_geneids
write.csv(lasso_geneids,file = 'Low_BMD/lasso.gene.csv',row.names = F)

beta_min <- as.matrix(coef(res.lasso, s = "lambda.min"))
beta_df <- data.frame(
  gene = rownames(beta_min),
  coef = as.numeric(beta_min)
)
beta_df <- beta_df[beta_df$coef != 0, ]
beta_df <- beta_df[-1,]
max_coef <- max(beta_df$coef) * 1.2

p <- ggplot(beta_df, aes(x = reorder(gene, coef), y = coef, fill = coef > 0)) +
  geom_col(width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#b5aa0f", "FALSE" = "#7ba39d")) +
  geom_label(aes(x = gene, y = max_coef, label = coef),
             hjust = 0, size = 3.5,
             fill = "#f7c173", color = "black", label.size = 0.3, label.padding = unit(0.15, "lines"))+
  
  labs(x = "Gene", y = "Coefficient", title = "OP Cohort: Gene Coefficient - GSE56814") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(min(beta_df$coef), max_coef * 2))
p

ggsave("Low_BMD/03.Coef.Bar.png", p, width = 8, height =5, dpi = 300,units = 'in',bg = 'white')
ggsave("Low_BMD/03.Coef.Bar.pdf", p, width = 8, height =5)
