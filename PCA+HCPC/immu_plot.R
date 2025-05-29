.libPaths("D:/PackagesR/R_LIBS")
rm(list=ls())
#运行完IOBR再运行
load(file = "clust_MR.Rdata")
library(ggplot2)
library(tidyverse)
library(corrplot)
library(pheatmap)
library(ggpubr)
library(reshape2)
library(ggsci)

colnames(im_xcell)<-sub('_xCell', '', colnames(im_xcell))  
df<-as.data.frame(im_xcell)
rownames(df)<-df[,1]
df<-df[,-1]
df<-as.data.frame(t(df))
df<- df[, Sam_order]
#df<-log2(df+1)

Type=annCol1
summary(factor(Type[,1]))
#cols<-colorRampPalette(c("red","white", "blue"))(95)
pdf("pheatmap_ssGSEA.pdf",width=12,height=8)
pheatmap(df,
         show_colnames = F,
         annotation_col = Type,
         #color = colorRampPalette(c("blue", "white","red"))(95),
         color = colorRampPalette(c( "white","red"))(95),
         cluster_cols = FALSE, 
         fontsize = 10)
dev.off()

library(tidyverse)
rt1=as.data.frame(t(df))
dat <- rt1 %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Complement_cell,value = Expression,-Sample)
dat$Cluster = Type$Cluster
library(ggrepel)
library(ggpubr)  
ggplot(dat,aes(Complement_cell,Expression,fill =Cluster)) + 
  geom_boxplot(outlier.shape = 21,color = "white") + 
  labs(title="Three Cluster of immune infiltration by xCell",x = "", y = "Enrichment scores") +
  theme(axis.text.x = element_text(angle=60,vjust = 0.5))+
  theme_classic()  + 
  theme(axis.title.y =  element_text(size = 10,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 12,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 10))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "right") + 
  scale_fill_manual(values = c("greenyellow","tomato2","deepskyblue1"))+
  stat_compare_means(aes(group=Cluster,label = ..p.signif..),method ="t.test",vjust =0.01,hide.ns = TRUE)


ggsave("boxplot_xCell.pdf",width =20,height =6)

genes<-c("USP8",'USP25',"USP21",'UFM1','UBL4A','UBE2V1','UBE2N',
         'UBE2L3','UBE2G2','UBE2B','OTUB2','ATG7','ISG15')
genes<-c("ISG15","UBE2B","UBE2V1")
expr_coad<-data.table::fread("GSE157103_tpm_COV.txt")
expr_coad<-as.data.frame(expr_coad)
rownames(expr_coad)<-expr_coad[,1]
expr_coad<-expr_coad[,-1]
genes_expr <- as.data.frame(t(expr_coad[rownames(expr_coad) %in% genes,]))
#

test<-im_ssgsea

genes_expr <- genes_expr[match(test$ID,rownames(genes_expr)),]
identical(test$ID,rownames(genes_expr))
library(linkET)
# install.packages("devtools")
#devtools::install_github("Hy4m/linkET", force = TRUE)
#packageVersion("linkET")
cor_res <- correlate(genes_expr, test[,-1],method = "spearman")
qcorrplot(cor_res) +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))
df_r <- cor_res$r %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(-1,names_to = "cell_type",values_to = "correlation")

df_p <- cor_res$p %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(-1,names_to = "cell_type",values_to = "pvalue")

df_cor <- df_r %>% 
  left_join(df_p) %>% 
  mutate(stars = cut(pvalue,breaks = c(-Inf,0.05,0.01,0.001,Inf),right = F,labels = c("***","**","*"," ")))
## Joining with `by = join_by(gene, cell_type)`

head(df_cor)
library(ggplot2)

ggplot(df_cor, aes(cell_type,gene))+
  geom_tile(aes(fill=correlation))+
  geom_text(aes(label=stars), color="black", size=4)+
  scale_fill_gradient2(low='#67B26F', high='#F2AA9D',mid = 'white',
                       limit=c(-1,1),name=paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_text(size=8,angle = 45,hjust = 1,color = "black"),
        axis.text.y = element_text(size=8,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank())
#https://aitechtogether.com/article/16734.html
#pheatmap画法
df_r <- cor_res$r %>% 
  as.data.frame() 
df_p <- cor_res$p %>% 
  as.data.frame() 

df_p[df_p<0.001]<-"***"
df_p[df_p>0.001 & df_p<0.01]<-"**"
df_p[df_p>0.01 & df_p<0.05]<-"*"
df_p[df_p>0.05]<-""



pheatmap(df_r, 
         show_colnames = TRUE,   # 是否显示列名
         show_rownames=TRUE,     # 是否显示行名
         fontsize=10,             # 字体大小
         color = colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(50), # 指定热图的颜色
         annotation_legend=TRUE, # 是否显示图例
         border_color=NA,        # 边框颜色 NA表示没有
         scale="none",           # 指定归一化的方式。"row"按行归一化，"column"按列归一化，"none"不处理
         cluster_rows = TRUE,    # 是否对行聚类
         cluster_cols = FALSE, 
         clustering_distance_rows = "correlation",
         cutree_rows =3,
         display_numbers = df_p,
         angle_col = 45,
         legend_breaks=c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8)
)

mat<-df_r
row_clust <- hclust(dist(mat,method = "euclidean"),method = "complete") # 行(特征)聚类
rowInd <- row_clust$order # 行(特征)的顺序

mat<- mat[rowInd,] # 将数据按照聚类结果重新排序
melt_mat <- melt(mat) # 融合数据,使之适用ggplot
library(ggtree)
h <- ggtree(row_clust,layout = "rectangular",branch.length = "none") # 行聚类树
p<-ggplot(df_cor, aes(cell_type,gene))+
  geom_tile(aes(fill=correlation))+
  geom_text(aes(label=stars), color="black", size=4)+
  scale_fill_gradient2(low='#67B26F', high='#F2AA9D',mid = 'white',
                       limit=c(-1,1),name=paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_text(size=8,angle = 45,hjust = 1,color = "black"),
        axis.text.y = element_text(size=8,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank())


p.all <- p%>% insert_left(h,width = 0.15) 
