#https://mp.weixin.qq.com/s/JqO7rVBMGGmOXRA8w8nDSg

.libPaths("D:/PackagesR/R_LIBS")
rm(list = ls())

# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

depens<-c('tibble', 'survival', 'survminer', 'sva', 'limma', "DESeq2","devtools",
          'limSolve', 'GSVA', 'e1071', 'preprocessCore', 'ggplot2', "biomaRt",
          'ggpubr', "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor", "timeROC","pracma")

for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(depen,update = FALSE)
}
#if (!requireNamespace("IOBR", quietly = TRUE))
 # devtools::install_github("IOBR/IOBR")

library(IOBR)

rt<-data.table::fread("GSE157103_tpm_COV.txt")
mrna_expr_tpm<-as.data.frame(rt)
rownames(mrna_expr_tpm)<-mrna_expr_tpm[,1]
mrna_expr_tpm<-mrna_expr_tpm[,-1]
expr_coad <- log2(mrna_expr_tpm+0.1)
expr_coad <- expr_coad[apply(expr_coad,1,sd)>0.5,]
dim(expr_coad)
tme_deconvolution_methods

# MCPcounter
im_mcpcounter <- deconvo_tme(eset = expr_coad,
                             method = "mcpcounter"
)
im_epic <- deconvo_tme(eset = expr_coad,
                       method = "epic",
                       arrays = F
)
# xCell
im_xcell <- deconvo_tme(eset = expr_coad,
                        method = "xcell",
                        arrays = F
)
# CIBERSORT
im_cibersort <- deconvo_tme(eset = expr_coad,
                            method = "cibersort",
                            arrays = F,
                            perm = 1000
)
save(im_cibersort,file = "GSE157103_COV_TMP_cobersort.Rdata")
load(file = "GSE157103_COV_TMP_cobersort.Rdata")
#perm：表示置换次数，数字越大运行时间越长，一般文章都设置为1000；
# QN：如果为芯片数据这里设为“T”
im_cibersort<-im_cibersort[,-c(16,24:26)]
# IPS
im_ips <- deconvo_tme(eset = expr_coad,
                      method = "ips",
                      plot = F
)
# quanTIseq
im_quantiseq <- deconvo_tme(eset = expr_coad,
                            method = "quantiseq",
                            scale_mrna = T
)
# ESTIMATE
im_estimate <- deconvo_tme(eset = expr_coad,
                           method = "estimate"
)
# TIMER
im_timer <- deconvo_tme(eset = expr_coad
                        ,method = "timer"
                        ,group_list = rep("coad",dim(expr_coad)[2])
)
# ssGSEA
load(file = "ssGSEA28.Rdata")
im_ssgsea <- calculate_sig_score(eset = expr_coad
                                 , signature = cellMarker # 这个28种细胞的文件需要自己准备
                                 , method = "ssgsea" # 选这个就好了
)
tme_combine <- im_mcpcounter %>% 
  inner_join(im_epic, by="ID") %>% 
  inner_join(im_xcell, by="ID") %>% 
  inner_join(im_cibersort, by="ID") %>% 
  inner_join(im_ips, by= "ID") %>% 
  inner_join(im_quantiseq, by="ID") %>% 
  inner_join(im_timer, by= "ID")%>% 
  inner_join(im_ssgsea, by= "ID")

# 绘制热图
# 查看各种结果的数值分布范围
range(im_epic[c(2, ncol(im_epic))]); range(im_xcell[c(2, ncol(im_xcell))]); range(im_cibersort[c(2, ncol(im_cibersort))])
# 数值分布范围较大，为了使热图更好看，需要先对数据进行处理
tme_combine1 <-  column_to_rownames(tme_combine, "ID")
range(tme_combine1)

# scale默认对列进行均一化，这里我们是要对不同算法结果，也就是行，进行均一化，因此需要转置
tme_combine2 <- scale(tme_combine1,
                      center = TRUE, scale = TRUE)

range(tme_combine2, na.rm = TRUE)

pheatmap(t(tme_combine2), 
         border_color = NA, 
         color = colorRampPalette(c("blue", "white","red"))(100), 
         cluster_rows = FALSE, # 行不聚类
         cluster_cols = TRUE, # 列不聚类
         show_rownames = FALSE, # 不显示行名
         show_colnames = FALSE)
# 还是存在极端值，因此对数据再次进行处理
tme_combine2[tme_combine2 > 5] = 5
tme_combine2[tme_combine2 < (-5)] = -5
pheatmap(t(tme_combine2), 
         border_color = NA, 
         color = colorRampPalette(c("blue", "white","red"))(100), 
         cluster_rows = FALSE, # 行不聚类
         cluster_cols = TRUE, # 列不聚类
         show_rownames = FALSE, # 不显示行名
         show_colnames = FALSE)
# 添加注释
Type <- rbind(data.frame(ID = colnames(im_ips)[2:ncol(im_ips)], method = "IPS"),
              data.frame(ID = colnames(im_timer)[2:ncol(im_timer)], method = "TIMER"),data.frame(ID = colnames(im_epic)[2:ncol(im_epic)], method = "EPIC"),
              data.frame(ID = colnames(im_mcpcounter)[2:ncol(im_mcpcounter)], method = "MCP_counter"),
              data.frame(ID = colnames(im_quantiseq)[2:ncol(im_quantiseq)], method = "Quanti-seq"),
              data.frame(ID = colnames(im_cibersort)[2:ncol(im_cibersort)], method = "CIBERSORT"),
              data.frame(ID = colnames(im_ssgsea)[2:ncol(im_ssgsea)], method = "ssGSEA"),
              data.frame(ID = colnames(im_xcell)[2:ncol(im_xcell)], method = "xCell"))
table(Type$method)
Type <- column_to_rownames(Type, "ID")

load(file = "clust_MR.Rdata")
for (i in 1:nrow(clust_mr)) {
  if(clust_mr[i,2]=="1")
    clust_mr[i,2]<-"Cluster 1"
  if(clust_mr[i,2]=="2")
    clust_mr[i,2]<-"Cluster 2"
  if(clust_mr[i,2]=="3")
    clust_mr[i,2]<-"Cluster 3"
}
colnames(clust_mr)<-c("ID","Cluster")
clust_mr<-as.data.frame(clust_mr)
annCol <- clust_mr[order(clust_mr$Cluster, decreasing = F), ]

rownames(annCol)<-annCol[,1]
annCol<-as.data.frame(annCol)

# 创建注释信息
# 列注释，位于热图顶端
annCol1<-annCol1 %>%
  arrange(Cluster)

# 行注释，位于热图左侧
annRow <- Type 

# 设置颜色
library(RColorBrewer)
display.brewer.all()
brewer.pal(7, "Set1")
annColors <- list(method = c("IPS"="pink",
                             "TIMER" = "#FFFF33",
                             "EPIC" = "#377EB8",
                             "MCP_counter" = "#E41A1C", 
                             "Quanti-seq" = "#FF7F00",
                             "CIBERSORT" = "#984EA3",
                             "ssGSEA" = "#A65628",
                             "xCell" = "#4DAF4A"),
                  Cluster= c("Cluster 1"="lightgreen","Cluster 2"="red","Cluster 3"="lightblue"))

# 按TCGA_NCPS排序


Sam_order <- rownames(annCol[order(annCol$Cluster),])
plot_data <- t(tme_combine2)[rownames(annRow), Sam_order]
annCol <- annCol[Sam_order, ,drop = F]
annCol1 <- annCol[,2,drop = F]

# pheatmap绘图
pdf("IOBR3.pdf", 10, 13)
pheatmap(plot_data,
         color = colorRampPalette(c("blue", "white","red"))(100), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         fontsize_row = 5,
         annotation_col = annCol1, 
         annotation_row = annRow,
         annotation_colors = annColors, 
         show_rownames = TRUE, # 显示行名
         show_colnames = FALSE,
         gaps_col = cumsum(table(annCol1$Cluster)), # 列分割
         gaps_row = c(6,12,20,30,41,62,90,157)# 行分割
)
dev.off()
a<- cumsum(table(annRow$method))
