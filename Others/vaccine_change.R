
rm(list=ls())  #清空环境内变量
.libPaths("D:/PackagesR/R_LIBS")

rt1<-data.table::fread("GSE220682_day0-22.csv")
rt2<-data.table::fread("GSE220682_day180-360.csv")
gene<-c("UBE2V1",  "UBE2L3" ,"UBE2B" ,"ISG15"  )
library(dplyr)

rt_gene1 <- rt1 %>% filter(GENE %in% gene)
rtlog1<-log(rt_gene1[,-1]+1)
rt_gene1[,c(2:ncol(rt_gene1))]<-rtlog1
write.table(rt_gene1,file = "GSE220682-day0-22_log.txt",sep = "\t",row.names = FALSE,
            col.names = TRUE)

rt_gene2 <- rt2 %>% filter(GENE %in% gene)
rtlog2<-log(rt_gene2[,-1]+1)
rt_gene2[,c(2:ncol(rt_gene2))]<-rtlog2
write.table(rt_gene2,file = "GSE220682-day180-360_log.txt",sep = "\t",row.names = FALSE,
            col.names = TRUE)
#筛选

rt3<-data.table::fread("GSE220682-day0-22_need.txt")
rt4<-data.table::fread("GSE220682-day180-360_need.txt")

rt5<-rbind(rt3,rt4)
order<-c("T1","T2","T3","T4")
dt <- rt5 %>%
  mutate(Stage = factor(Stage, levels = order)) %>%
  arrange(Stage)

library("ggpubr")
library("ggsci")
library("scales")
dir="./plot2/"  
#ISG15
#UBE2B
#UBE2L3
#UBE2V1
##########boxplot4##########
pal=pal_simpsons("springfield", 
                 alpha = 0.2)(4) 
pal;show_col(pal)                     #查看颜色

#改
outFile=paste0(dir,"UBE2L3.pdf" )   #输出文件

#设置分组
x=colnames(dt)[2]
#改
y=colnames(dt)[5]
Stage=levels(factor(dt$Stage))
dt$Stage=factor(dt$Stage, levels=Stage)
comp=combn(Stage,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
#绘制
pdf(file=outFile, width=5.5, height=5)
#改
boxplot=ggboxplot(dt, x="Stage", y="UBE2L3", color="Stage",
                  xlab=x,
                  ylab=y,
                  legend.title=x,
                  add = "jitter")+ 
  scale_colour_manual(values = c("#FED43999","#709AE199","#8A919799","#D2AF8199"))+
  stat_compare_means(comparisons = my_comparisons)
print(boxplot)
dev.off()


