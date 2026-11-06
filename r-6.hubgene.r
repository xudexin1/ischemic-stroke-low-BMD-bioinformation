rm(list = ls())
setwd('/home/data/project/65.CCPJ-2003-6/')
file.dir = '06_hubgene'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)
lasso_is <- read.csv('../05_Lasso/IS/lasso.gene.csv')
lasso_op <- read.csv('../05_Lasso/Low_BMD/lasso.gene.csv')
hubgene <- data.frame(symbol = intersect(lasso_is$x,lasso_op$x))
write.csv(hubgene,file = 'hubgene.csv',row.names = F,quote = F)


venn_list <- list('Lasso - IS' = lasso_is$x, 'Lasso - low BMD' = lasso_op$x)
text <- paste0("**HubGenes:**<br>", paste(hubgene$symbol, collapse = "<br>"))
venn <- ggvenn(venn_list, 
               fill_color = c("#BB5566", "#0072B2"), 
               stroke_size = 0.5, set_name_size = 5) +
  theme(plot.margin = margin(10, 50, 10, 10))

venn
box <- ggplot() +
  geom_textbox(
    aes(x = 0.5, y = 0.5, label = text),  # 中心位置
    box.color = "#BB5566",
    fill = "white",
    width = unit(0.7, "npc"),   # 宽度占满整个画布
    height = unit(0.7, "npc"),  # 高度占满整个画布（新加）
    hjust = 0.5, vjust = 0.5, # 居中定位文本框
    halign = 0,               # 文本框内靠左对齐文字
    size = 5,
    lineheight = 1.5
  ) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))
box
# patchwork组合
final_plot <- (venn | box) +
  plot_layout(widths = c(2, 1), heights = c(1, 1)) &
  theme(plot.margin = margin(5, 5, 5, 5))
final_plot
# 保存图像，防止裁切
ggsave("01.venn_plot.pdf", plot = final_plot,
       width = 8, height = 4)
ggsave("01.venn_plot.png", plot = final_plot,
       width = 8, height = 4, units = "in", dpi = 300)

