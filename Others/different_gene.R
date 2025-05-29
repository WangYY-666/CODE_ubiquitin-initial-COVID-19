
rm(list=ls())  #清空环境内变量
.libPaths("D:/PackagesR/R_LIBS")

#rt1<-data.table::fread("GSE161731_counts.csv")
rt2<-data.table::fread("GSE161731_clinical_data.txt")
rt3<-data.table::fread("GSE161731_xpr_tpm_geo.txt")

gene<-c("UBE2V1",  "UBE2L3" ,"UBE2B" ,"ISG15"  )
library(dplyr)

rt_gene <- rt3 %>% filter(Gene %in% gene)
rtlog<-log(rt_gene[,-1]+1)
rt_gene[,c(2:ncol(rt_gene))]<-rtlog
rt4<-as.data.frame(t(rt_gene))
rt4<-cbind(id=rownames(rt4),rt4)
colnames(rt4)<-rt4[1,]
rt4<-rt4[-1,]
colnames(rt4)[1]<-"id"
rt<-merge(rt4,rt2,by="id")

dt1<-rt[,c(1,9)]
dt1[c(1:19),2]<-"healthy"
colnames(dt1)[2]<-"group"
order <- c("healthy", "early","middle" , "late")

library("ggpubr")
library("ggsci")
library("scales")
dir="./plot/"  
#ISG15
dt_1<-merge(dt1,rt[,c(1,2)],by="id")
dt_1$ISG15<-as.numeric(dt_1$ISG15)

#UBE2B
dt_1<-merge(dt1,rt[,c(1,3)],by="id")
dt_1$UBE2B<-as.numeric(dt_1$UBE2B)

#UBE2V1
dt_1<-merge(dt1,rt[,c(1,4)],by="id")
dt_1$UBE2V1<-as.numeric(dt_1$UBE2V1)

#UBE2L3
dt_1<-merge(dt1,rt[,c(1,5)],by="id")
dt_1$UBE2L3<-as.numeric(dt_1$UBE2L3)

dt <- dt_1 %>%
  mutate(group = factor(group, levels = order)) %>%
  arrange(group)


##########boxplot4##########
 
#改
outFile=paste0(dir,"UBE2L3.pdf" )   #输出文件
pal=pal_simpsons("springfield", 
                 alpha = 0.6)(4) 

pal;show_col(pal)                     #查看颜色

#设置分组
x=colnames(dt)[2]
y=colnames(dt)[3]
group=levels(factor(dt$group))
dt$group=factor(dt$group, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制
pdf(file=outFile, width=5.5, height=5)
#改
boxplot=ggboxplot(dt, x="group", y="UBE2L3", color="group",
                  xlab=x,
                  ylab=y,
                  legend.title=x,
                  add = "jitter")+ 
  scale_colour_manual(values = c("#FED43999","#709AE199","#8A919799","#D2AF8199"))+
  stat_compare_means(comparisons = my_comparisons)
print(boxplot)
dev.off()


