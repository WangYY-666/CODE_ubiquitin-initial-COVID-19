.libPaths("D:/PackagesR/R_LIBS")
rm(list = ls())
library(limma)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(ggplot2)
library(RColorBrewer)

load(file = "GSE157103_C_normalized.Rdata")
#c("ISG15","UBE2B","UBE2L3","UBE2V1")
gene<-"UBE2V1"
cutoff<-rowMeans(rt3[gene,])
rt<-rt3
group <- ifelse(rt[gene,]>cutoff, "High", "Low")   
group <- factor(group,levels = c("High","Low"))
table(group)
#差异分析
#所选基因高表达vs低表达
#logfc为正说明（eg.ISG15表达高时该基因表达高）
design <- model.matrix(~group)
colnames(design) <- levels(group)
fit <- lmFit(rt,design)
fit2 <- eBayes(fit)
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf)

importgene<-allDiff[allDiff$P.Value<0.05 & abs(allDiff$logFC)>1,]
geneup<-allDiff[allDiff$P.Value<0.05 & allDiff$logFC>1,]
genedown<-allDiff[allDiff$P.Value<0.05 & allDiff$logFC < -1,]
#差异基因集提取：
up <- rownames(geneup)#差异上调
down <- rownames(genedown)#差异下调
diff <- c(up, down)#所有差异基因
#ID转换
#查看可转换的ID类型：
columns(org.Hs.eg.db)
#以总差异基因集为例进行转换：
##上调或下调基因的单独转换相同，这里不再赘述
diff_entrez <- bitr(diff,
                    fromType = "SYMBOL",#现有的ID类型
                    toType = "ENTREZID",#需转换的ID类型
                    OrgDb = "org.Hs.eg.db")
head(diff_entrez) #提示少量无法映射，不同数据库间ID转换存在少量缺失是正常现象

#KEGG富集分析(超几何分布检验)：
KEGG_diff <- enrichKEGG(gene = diff_entrez$ENTREZID,
                        organism = "hsa", #物种Homo sapiens (智人)
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        minGSSize = 10,
                        maxGSSize = 500)

#将ENTREZ重转为symbol：
KEGG_diff <- setReadable(KEGG_diff,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID")
#View(KEGG_diff@result)
#colnames(KEGG_diff@result) #解释不同列名含义
#计算Rich Factor（富集因子）：
KEGG_diff2 <- mutate(KEGG_diff,
                     RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

#计算Fold Enrichment（富集倍数）：
KEGG_diff2 <- mutate(KEGG_diff2, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
KEGG_diff2@result$RichFactor[1:6]
KEGG_diff2@result$FoldEnrichment[1:6]
#提取KEGG富集分析结果表：
KEGG_result <- KEGG_diff2@result

#保存富集结果到本地：
#save(KEGG_diff2, KEGG_result, file = c("KEGG_diff.Rdata"))
write.csv(KEGG_result, file = c('KEGG_UBE2V1_diff.csv'))
#查看感兴趣KEGG pathway：
#browseKEGG(KEGG_diff, "hsa04610") #直接跳转到KEGG
#富集条形图绘制：
