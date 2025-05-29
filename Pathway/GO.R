.libPaths("D:/PackagesR/R_LIBS")
rm(list = ls())
library(limma)

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


#分开分析
geneup<-allDiff[allDiff$P.Value<0.05 & allDiff$logFC>1,]
genedown<-allDiff[allDiff$P.Value<0.05 & allDiff$logFC < -1,]


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

rt<-rownames(genedown)
genes<-as.vector(rt)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)       
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
colnames(out)=c("symbol","entrezID")
out=out[is.na(out[,"entrezID"])==F,]
out<-as.data.frame(out)
gene=out$entrezID

ego_ALL <- enrichGO(gene = gene,
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)

ego_CC <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_BP <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

ego_MF <- enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)


ego_ALL <- as.data.frame(ego_ALL)
ego_result_BP <- as.data.frame(ego_BP)
ego_result_CC <- as.data.frame(ego_CC)
ego_result_MF <- as.data.frame(ego_MF)

ego_result_BP <- ego_result_BP[order(ego_result_BP$Count,decreasing = T),]
ego_result_CC <- ego_result_CC[order(ego_result_CC$Count,decreasing = T),]
ego_result_MF <- ego_result_MF[order(ego_result_MF$Count,decreasing = T),]
ego <- rbind(ego_result_BP,ego_result_CC,ego_result_MF)
write.csv(ego_ALL,file = "UBE2V1/ego_ALL.csv",row.names = T)
write.csv(ego_result_BP,file = "UBE2V1/ego_result_BP.csv",row.names = T)
write.csv(ego_result_CC,file = "UBE2V1/ego_result_CC.csv",row.names = T)
write.csv(ego_result_MF,file = "UBE2V1/ego_result_MF.csv",row.names = T)
write.csv(ego,file = "UBE2V1/ego.csv",row.names = T)



