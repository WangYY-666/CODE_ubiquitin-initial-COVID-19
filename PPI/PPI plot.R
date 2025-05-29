
rm(list=ls())
####----load R Package----####
library(tidyverse)
library(Biostrings)
library(ggfun)
library(ggraph)
library(tidygraph)
library(purrr)
library(rlang)
library(patchwork)
source("get_node_edge.R")
####----load Data----####
# STRING output
ppi_df <- read_delim(file = "0.4/string_interactions_short.tsv", delim = "\t", col_names = T)


edge_out <- ppi_df[,c(1,2,13)]%>%rename_with(~c("from","to","Combined_score"),c(1,2,3))
node_out <- ppi_df[,c(1,2)]%>%rename_with(~c("from","to"),c(1,2))

graph_df <- as_tbl_graph(edge_out) %>%
  tidygraph::mutate(Popularity = centrality_degree(mode = 'out'))


p1 <- ggraph(graph_df, layout = 'linear', circular = TRUE) + 
  geom_edge_arc(aes(color = Combined_score), alpha = 1) + 
  geom_node_point(aes(size = Popularity), fill= "#7fcdbb",color = "#000000",alpha=1, shape = 21) + 
  # 新增节点标签层
  geom_node_text(
    aes(label = name),          # 使用节点的 name 列作为标签
    repel = TRUE,                # 自动避开重叠
    size = 3,                   # 字体大小
    color = "black"            # 字体颜色        
  ) + 
  scale_edge_colour_distiller(palette = "RdPu", direction = 1) + 
  coord_fixed() +
  scale_size(range = c(5,15)) +
  theme_nothing()

p1

library(ggraph)
library(tidygraph)

p1 <- ggraph(graph_df, layout = 'linear', circular = TRUE) + 
  geom_edge_arc(aes(color = Combined_score), 
                alpha = 0.6, 
                edge_width = 0.8) +
  scale_edge_color_distiller(
    palette = "RdPu", 
    direction = 1,
    name = "Interaction Score"
  ) +
  
  # 修改点：添加 Popularity 的图例
  geom_node_point(
    aes(size = Popularity),
    fill = "#4D96B7",
    color = "#2D4059",
    alpha = 0.9,
    shape = 21,
    stroke = 0.6
  ) +
  scale_size_continuous(
    range = c(4, 12),
    name = "Node Degree (Popularity)",  # 图例标题
    breaks = c(min(graph_df %N>% pull(Popularity)), 
               median(graph_df %N>% pull(Popularity)), 
               max(graph_df %N>% pull(Popularity))),  # 自定义刻度点
    labels = c("Low", "Medium", "High")  # 图例标签
  ) +
  
  geom_node_text(
    aes(label = name),
    repel = TRUE, 
    size = 3.5, 
    color = "#2D4059",
    family = "sans",
    box.padding = 1.5,
    point.padding = 0.5,
    segment.color = "grey70",
    segment.size = 0.3,
    max.overlaps = 20
  ) +
  
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.box = "horizontal",  # 水平排列两个图例
    legend.spacing.x = unit(0.5, "cm")  # 图例间距
  )

print(p1)
