.libPaths("D:/PackagesR/R_LIBS")
rm(list = ls())

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

#ISG15 UBE2B UBE2L3 UBE2V1
ego_result_BP <- data.table::fread("UBE2V1/GO/ego_result_BP.csv")
ego_result_CC <- data.table::fread("UBE2V1/GO/ego_result_CC.csv")
ego_result_MF <- data.table::fread("UBE2V1/GO/ego_result_MF.csv")

ego_result_BP <- ego_result_BP[order(ego_result_BP$Count,decreasing = T),]
ego_result_CC <- ego_result_CC[order(ego_result_CC$Count,decreasing = T),]
ego_result_MF <- ego_result_MF[order(ego_result_MF$Count,decreasing = T),]

ego_result_BP$Description
ego_result_CC$Description
ego_result_MF$Description

#ISG15
#ego_result_BP <- ego_result_BP[c(1:6,11:13,15,16,20,25), ]
#ego_result_CC <- ego_result_CC[1, ]
#ego_result_MF <- ego_result_MF[c(3,4,10), ]
#UBE2B
#ego_result_BP <- ego_result_BP[c(1:11,13), ]
#ego_result_CC <- ego_result_CC
#ego_result_MF <- ego_result_MF[c(1:4,9,12), ]
#UBE2L3
#ego_result_BP <- ego_result_BP[c(2,3,5,6,7,10,14:17,25,27,30,33,40), ]
#ego_result_CC <- ego_result_CC[4,]
#ego_result_MF <- ego_result_MF[c(2,9,12), ]
#UBE2V1
ego_result_BP <- ego_result_BP[c(3,5,7,9,13:17,19,20,39,40,52,62), ]
ego_result_CC <- ego_result_CC[c(10,15,16),]
ego_result_MF <- ego_result_MF[c(1,6,12,19,28), ]

go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),                         
  Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", nrow(ego_result_BP)), 
                rep("cellular component", nrow(ego_result_CC)),
                rep("molecular function", nrow(ego_result_MF))), 
              levels=c("biological process", "cellular component","molecular function" )))
for(i in 1:nrow(go_enrich_df)){
  description_splite=strsplit(go_enrich_df$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:5],collapse = " ") 
  go_enrich_df$Description[i]=description_collapse
  go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
}


go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
library(RColorBrewer)
COLS <- brewer.pal(name="Dark2",3)

ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
  geom_bar(stat="identity", width=0.9) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = "The GO Terms of UBE2V1")+
  theme_bw()
ggsave("UBE2V1/GO整体分析/GO.pdf",width =8,height =4)


