.libPaths("D:/PackagesR/R_LIBS")
rm(list = ls())

#引用需要用到的包，没有的话请提前进行安装
library(GSEABase)
library(GSVA)
library(limma)
library(pheatmap)
library(ggplot2)

#文件的读入与处理
#读入表达矩阵并处理
load(file = "GSE157103_C_normalized.Rdata")
rt=as.matrix(rt3)
exp=rt
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
rt=normalizeBetweenArrays(as.matrix(rt))


#读入工作目录下的c2.cp.kegg.v7.1.symbols.gmt文件
#keggSet <- getGmt("c2.cp.kegg.v7.4.symbols.gmt")
C7set<-getGmt("c7.all.v2023.2.Hs.symbols.gmt")
#开始进行GSVA分析
#kegg <- gsva(expr=as.matrix(rt), keggSet, kcdf="Gaussian",method = "gsva",parallel.sz=1)
#exprSet <- kegg
#write.csv(kegg,"kegg__C7_results.csv")
#save(kegg,file = "kegg.Rdata")
C7 <- gsva(expr=as.matrix(rt), C7set, kcdf="Gaussian",method = "gsva",parallel.sz=1)
exprSet <- C7
write.csv(C7,"gsva__C7_results.csv")
save(C7,file = "C7.Rdata")

#load(file = "kegg.Rdata")
#("ISG15","UBE2B","UBE2L3","UBE2V1")
load(file = "C7.Rdata")
gene<-"UBE2V1"
cutoff<-rowMeans(rt3[gene,])
group <- ifelse(rt[gene,]>cutoff, "High", "Low")   
group <- factor(group,levels = c("High","Low"))
table(group)
#差异分析
design <- model.matrix(~group)
colnames(design) <- levels(group)
fit <- lmFit(exprSet,design)
fit2 <- eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf)
allDiff1<-rbind(pathway=colnames(allDiff),allDiff)
#write.table(allDiff1,file="UBE2V1_KEGG_alldiff.xls",sep="\t",quote=F,col.names = F)
write.table(allDiff1,file="UBE2V1_C7_alldiff.xls",sep="\t",quote=F,col.names = F)

orderdiff<-allDiff[order(abs(allDiff$logFC),decreasing = T),]
index <- rownames(orderdiff)[orderdiff$adj.P.Val<0.05 ]
annotation_col <- data.frame(group)

