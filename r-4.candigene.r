rm(list = ls())
setwd('/home/data/project/65.CCPJ-2003-6/')
file.dir = '04_candigene'
if (! dir.exists(file.dir)){
  dir.create(file.dir)
}
setwd(file.dir)

deg1 <- read.csv('../01_DEGs_IS/DEG_sig.csv',row.names = 1)
deg2 <- read.csv('../02_DEGs_OP/DEG_sig.csv',row.names = 1)
up1 <- subset(deg1,change == 'UP')
up2 <- subset(deg2,change == 'UP')
down1 <- subset(deg1,change == 'DOWN')
down2 <- subset(deg2,change == 'DOWN')

up <- data.frame(symbol = intersect(rownames(up1),rownames(up2)),change = 'Up')
down <- data.frame(symbol = intersect(rownames(down1),rownames(down2)),change = 'Down')
common <- rbind(up,down)

write.csv(common,file = 'candigene.csv',row.names = F,quote = F)

library(ggvenn)
library(ggtext)
library(patchwork)
up_list <- list('IS Up' = rownames(up1), 'Low BMD Up' = rownames(up2))
down_list <- list('IS Down' = rownames(down1), 'Low BMD Down' = rownames(down2))

up_common <- intersect(up_list$`IS Up`, up_list$`Low BMD Up`)
down_common <- intersect(down_list$`IS Down`, down_list$`Low BMD Down`)

up_text <- paste0("**Up-regulated genes:**<br>", paste(up_common, collapse = "<br>"))
down_text <- paste0("**Down-regulated genes:**<br>", paste(down_common, collapse = "<br>"))
venn_up <- ggvenn(up_list, 
                  fill_color = c("#56B4E9", "#009E73"), 
                  stroke_size = 0.5, set_name_size = 5) +
  theme(plot.margin = margin(10, 50, 10, 10))

venn_down <- ggvenn(down_list, 
                    fill_color = c("#F0E442", "#D55E00"), 
                    stroke_size = 0.5, set_name_size = 5) +
  theme(plot.margin = margin(10, 50, 10, 10))

up_box <- ggplot() +
  geom_textbox(
    aes(x = 0.5, y = 0.5, label = up_text),  # 中心位置
    box.color = "#009E73",
    fill = "white",
    width = unit(0.7, "npc"),   # 宽度占满整个画布
    height = unit(0.7, "npc"),  # 高度占满整个画布（新加）
    hjust = 0.5, vjust = 0.5, # 居中定位文本框
    halign = 0,               # 文本框内靠左对齐文字
    size = 4,
    lineheight = 1.1
  ) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))
up_box
down_box <- ggplot() +
  geom_textbox(
    aes(x = 0.5, y = 0.5, label = down_text),  # 中心位置
    box.color = "#D55E00",
    fill = "white",
    width = unit(0.7, "npc"),   # 宽度占满整个画布
    height = unit(0.7, "npc"),  # 高度占满整个画布（新加）
    hjust = 0.5, vjust = 0.5, # 居中定位文本框
    halign = 0,               # 文本框内靠左对齐文字
    size = 4,
    lineheight = 1.1
  ) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))
down_box

# patchwork组合
final_plot <- (venn_up | up_box) / (venn_down | down_box) +
  plot_layout(widths = c(2, 1), heights = c(1, 1)) &
  theme(plot.margin = margin(5, 5, 5, 5))
final_plot
# 保存图像，防止裁切
ggsave("01.venn_plot.pdf", plot = final_plot,
       width = 10, height = 8)
ggsave("01.venn_plot.png", plot = final_plot,
       width = 10, height = 8, units = "in", dpi = 300)
