
.libPaths("D:/PackagesR/R_LIBS")
rm(list=ls())

library(factoextra) # 可视化
library(FactoMineR) # PCA 聚类
library(tidyverse)  # %>%
library(ggcorrplot)
library(tidyverse)  
library(pca3d)

rt<-data.table::fread("GSE157103_ec.txt")
rt<-as.data.frame(rt)
rownames(rt)<-rt[,1]
rt1<-rt[,-1]
rt1<-rt1[-1,]
rt2<-as.data.frame(t(rt1))
rt3<-as.data.frame(lapply(rt2,as.numeric))
rownames(rt3)<-rownames(rt2)
rt4<-rt3[-c(101:126),]
gene<-c("USP8",'USP25',"USP21",'UFM1','UBL4A','UBE2V1','UBE2N',
       'UBE2L3','UBE2G2','UBE2B','OTUB2','ATG7','ISG15')
rt5<-rt4[,gene]

#数据标准化
df <- rt5 %>%
  filter(rowSums(.) > 1) %>%
  na.omit() %>%
  mutate_all(~log2(. + 1)) %>%
  scale() %>%
  #t() %>%
  scale(center = TRUE, scale = FALSE)
dim(df)

# Compute PCA with ncp = 3
res.pca <- PCA(df, ncp = 3, graph = FALSE)

fviz_eig(res.pca, addlabels = TRUE)
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 13)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 13)
# Contributions of variables to PC3
fviz_contrib(res.pca, choice = "var", axes = 3, top = 13)
#all
fviz_contrib(res.pca, choice = "var", axes = c(1:3), top = 13)
fviz_pca_var(res.pca, 
             col.var="contrib",   # 根据贡献度着色       
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),          
             repel = TRUE)      
fviz_pca_ind(res.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
library(corrplot) 
corrplot(res.pca$var$contrib, is.corr=FALSE)
# Create a random continuous variable of length 10 
set.seed(123) 
my.cont.var <- rnorm(13) 
# Color variables by the continuous variable 
fviz_pca_var(res.pca, col.var = my.cont.var,              
             gradient.cols = c("blue", "yellow", "red"),              
             legend.title = "Cont.Var")
# Create a grouping variable using kmeans 
# Create 3 groups of variables (centers = 3) 
set.seed(123) 
res.km <- kmeans(res.pca$var$coord, centers = 3, nstart = 25) 
grp <- as.factor(res.km$cluster) 
# Color variables by groups 
fviz_pca_var(res.pca, col.var = grp,               
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),              
             legend.title = "Cluster")
#识别与给定主成分最显着相关的变量
res.desc <- dimdesc(res.pca, axes = c(1,2,3), proba = 0.05) # Description of dimension 1 res.desc$Dim.1
res.desc$Dim.1
res.desc$Dim.2
res.desc$Dim.3

fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color 
                col.ind = "#696969" # Individuals color 
)

# Compute hierarchical clustering on principal components
res.hcpc <- HCPC(res.pca, graph = FALSE)
#聚类树可视化
fviz_dend(res.hcpc,
          cex = 0.7,                     # Label size
          palette = "lancet",            # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "lancet",        # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)

#pca
clust<-as.data.frame(res.hcpc$data.clust)
clust_mr<-cbind(id=rownames(clust),clust=clust$clust)
save(clust_mr,file = "clust_MR.Rdata")
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = clust$clust, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Clusts"
)
fviz_cluster(res.hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE,  # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)
# Principal components + tree
plot(res.hcpc, choice = "3D.map")

res.hcpc$desc.axes$call$X
res.hcpc$desc.axes$call
res.hcpc$desc.axes$quanti
res.hcpc$desc.ind$para
res.hcpc$desc.axes$quanti.var
res.hcpc$desc.var$quanti.var

#每个聚类最多的定量变量
res.hcpc$desc.var$quanti



pca <- prcomp(df, scale. = TRUE)
get_eig(pca)
clust<-as.data.frame(res.hcpc$data.clust)
pca3d(pca, group = clust$clust, fancy = TRUE, 
      bg = "black", 
      show.group.labels = TRUE,
      axes.color = "white", new = TRUE)
library(rgl)  
writeWebGL(file = "my_rgl_plot.html")

pca3d(pca, group = clust$clust,  bg = "white", 
             show.ellipses = TRUE, ellipse.ci = 0.95, 
             show.plane = TRUE,show.scale = TRUE,
             show.shadows = TRUE)
rgl.snapshot(filename = "plot.png")

