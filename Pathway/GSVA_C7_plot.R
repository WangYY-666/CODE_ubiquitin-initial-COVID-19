.libPaths("D:/PackagesR/R_LIBS")
rm(list = ls())

library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(dplyr)

#ISG15 UBE2B UBE2L3 UBE2V1
df <- data.table::fread("ISG15/ISG15_C7_alldiff.xls")
df<-as.data.frame(df)
rownames(df)<-df[,1]
df<-df[,-1]
df$gene <- rownames(df)

cutoff<-0.3
df$threshold = factor(ifelse(df$P.Value  < 0.05 & abs(df$logFC) >= cutoff, ifelse(df$logFC >= cutoff ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
df$gene <- row.names(df)
df<-df[order(abs(df$logFC),decreasing = T),]
p1=ggplot(df,aes(x=logFC,y= -log10(P.Value),size = abs(logFC),fill = threshold))+
  geom_point(colour = "black", shape = 21, stroke = 0.5)+
  scale_size(limits  = c(0, 1))+
  scale_fill_manual(values=c("#fe0000","#13fc00","#bdbdbd"))+
  
  ylab('-log10 (Pvalue)')+
  xlab('log2 (FoldChange)')+
  geom_vline(xintercept=c(-cutoff,cutoff),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.5) +
  theme_classic( 
    base_line_size = 1 
  )+  
  
  guides(fill=guide_legend(override.aes = list(size=5)))+
  theme(axis.title.x = element_text(size = 10, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 10,
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.text = element_text(color="black", 
                                   size = 7, 
                                   face = "bold")
  ) 
for_label <- df %>% 
  filter(abs(logFC)>cutoff &  adj.P.Val< 0.05)

rownames(for_label)
#ISG15
#for_label<-for_label[c(32,34,46,48,49,53,56,58,73,76,97,101,104,107,118,123,124),]
#UBE2B
#for_label<-for_label[c(1,7,9,11:13,15,16,27:30,35,36,39,44,51,53,56,57,62,63,68,73,75:79,83,101,117),]
#for_label<-for_label[-c(14,15,18,19,24,28,32),]
#UBE2L3
#for_label<-for_label[c(16,17,26,27,30,34,41:43,49,51,64,66,78,80,83,85:88,106,111,115,121,123,124,129,131,135,151,166,170,171,178,179,181,183,188,194,200,202,204),]
#for_label<-for_label[-c(16,25,27,28,41),]
#UBE2V1
#for_label<-for_label[c(5,6,13,15,18,20:23,26,27,31,31,35,38,43,45,51,52,54,56,57,60,72,73,78:80,82,83,85:87,94,96,98,103,104,107,111,114,119,122,127,130,136,139,148,160:162,168,172,178,187,190,195:197,200,217,218,221,223,227:229,234,235,238),]
#for_label<-for_label[-c(11,13,14,16,18,20,31,33,34,38,40,41,43,45,56,63),]

library(stringr)
for_label$gene<-sub("^[^_]+_", "",for_label$gene)

#save(for_label,file = "ISG15/GSVA_C7_label_pathway.Rdata")
p1 +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = gene),
    data = for_label,alpha = 0.6,fontface="bold", 
    color="black"
  )+ theme_classic(base_size = 8)+
  theme(axis.title.y =  element_text(size = 20,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 20,angle = 45,vjust = 0.5,hjust = 0.5))+
  theme(axis.title.x =  element_text(size = 20,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.y = element_text(size = 20,angle = 45,vjust = 0.5,hjust = 0.5))+
  
  theme(legend.position = "NA")
ggsave("UBE2B/GSVA_C7.pdf", height = 12, width = 15)
