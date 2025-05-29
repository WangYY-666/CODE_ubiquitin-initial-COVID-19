.libPaths("D:/PackagesR/R_LIBS")
rm(list = ls())

rt<-data.table::fread("UBE2V1/KEGG/KEGG_UBE2V1_diff.csv")

#富集条形图绘制：
library(enrichplot)
library(ggplot2)
#以富集结果表Top20为例：
rt1<- rt[order(rt$pvalue,decreasing = F),]
rt1$Description

#ISG15
#rt2<-rt1[c(3,5,7:11,13:15,16,25,57),]
#UBE2B
#rt2<-rt1[c(1,3,4,6,7,11,12,31,85,112),]
#UBE2L3
#rt2<-rt1[c(1:4,12,14,15,17,18,20,26,27,29,61),]
#UBE2V1
rt2<-rt1[c(1,2,4,8,10,12,18,22,51),]
#指定绘图顺序（转换为因子）：
rt2$pathway <- factor(rt2$Description, levels = rev(rt2$Description))


#Top20富集数目条形图：
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11)) #自定义主题
#Top20显著富集条形图：
p1 <- ggplot(data = rt2,
             aes(x = -log10(pvalue), y = pathway, fill = Count)) +
  geom_bar(stat = "identity",width = 0.8) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(x = "-log10(pvalue)",
       y = "pathway",
       title = "KEGG enrichment barplot") +
  theme_bw() +
  mytheme
#x11()
p1
