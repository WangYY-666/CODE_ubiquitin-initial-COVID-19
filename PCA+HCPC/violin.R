
.libPaths("D:/PackagesR/R_LIBS")
library(ggpubr)
rm(list = ls())
rt2=read.table("violin.txt",sep="\t",header=T,check.names=F) 
load(file = "GSE157103_Clinic_tidy.Rdata")
load(file = "clust_MR.Rdata")
a_cli<-as.data.frame(clust_mr)
for (i in 1:nrow(a_cli)) {
  if(a_cli[i,2]=="1")
    a_cli[i,2]<-"Cluster 1"
  if(a_cli[i,2]=="2")
    a_cli[i,2]<-"Cluster 2"
  if(a_cli[i,2]=="3")
    a_cli[i,2]<-"Cluster 3"
}
rt2<-rtt
colnames(rt2)[1]<-"C"
colnames(a_cli)[1]<-"C"
colnames(a_cli)[2]<-"cluster"
rt<-merge(a_cli,rt2,by="C")
rt2<-rt
COLS <- c("lightgreen","tomato3", "deepskyblue3")

#APACHE Ⅱ评分最初设计为入ICU 24小时内最差值评分，因而一般不用于连续动态评价患者的病情危重程度。
rt3<-rt2[,c(3,6,7,9:18)]
rt3<-as.data.frame(lapply(rt3,as.numeric))
rt2<-cbind(rt2[,c(1,2)],rt3)
head(rt2)
dp <- ggplot(rt2, aes(x=cluster, y=HFDP45, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of HFD45",x="cluster", y = "HFD45")+ 
  stat_compare_means() # Add global p-value
p1=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.line=element_line(colour="black"))+
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")


head(rt2)
#rt2$ventilator_free_days<-log(rt2$ventilator_free_days+1)
dp <- ggplot(rt2, aes(x=cluster, y=ventilator_free_days, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of ventilator_free_days",x="cluster", y = "ventilator_free_days")+ 
  #labs(title="Plot of ventilator_free_days",x="cluster", y = "log(ventilator_free_days+1)")+ 
  stat_compare_means() 
p2=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")





head(rt2)
dp <- ggplot(rt2, aes(x=cluster, y=`ferritin.ng.ml.`, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of Ferritin",x="cluster", y = "Ferritin(ng/ml)")+ 
  stat_compare_means() 
p3=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")



head(rt2)
dp <- ggplot(rt2, aes(x=cluster, y=`crp..mg.l.`, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of CRP",x="cluster", y = "CRP (mg/l)")+ 
  stat_compare_means() 
p4=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")


head(rt2)
# 对value进行对数转换以减少离群点的影响  
rt2$ddimer..mg.l_feu.<-log(rt2$ddimer..mg.l_feu. + 1) 
dp <- ggplot(rt2, aes(x=cluster, y=`ddimer..mg.l_feu.`, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  #labs(title="Plot of ddimer",x="cluster", y = "ddimer (mg/l_feu)")+ 
  labs(title="Plot of ddimer",x="cluster", y = "log(ddimer+1)  (mg/l_feu)")+ 
  stat_compare_means()
p5=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")


head(rt2)
rt2$procalcitonin..ng.ml.<-log(rt2$procalcitonin..ng.ml.+1)
dp <- ggplot(rt2, aes(x=cluster, y=`procalcitonin..ng.ml.`, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  #labs(title="Plot of Procalcitonin",x="cluster", y = "Procalcitonin (ng/ml)")+ 
  labs(title="Plot of Procalcitonin",x="cluster", y = "log(Procalcitonin+1) (ng/ml)")+ 
  stat_compare_means()
p6=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")

dp <- ggplot(rt2, aes(x=cluster, y=`charlson.score`, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of charlson score",x="cluster", y = "score")+ 
  stat_compare_means() 
p7=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")

dp <- ggplot(rt2, aes(x=cluster, y=sofa, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of SOFA score",x="cluster", y = "score")+ 
  stat_compare_means() 
p8=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")

dp <- ggplot(rt2, aes(x=cluster, y=`lactate..mmol.l.`, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of lactate",x="cluster", y = "lactate (mmol/l)")+ 
  stat_compare_means() 
p9=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")

dp <- ggplot(rt2, aes(x=cluster, y=fibrinogen, fill=cluster)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of fibrinogen",x="cluster", y = "fibrinogen")+ 
  stat_compare_means() 
p10=dp + scale_fill_manual(values=COLS )+ theme_minimal()+ 
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")

library(cowplot)

plot_grid(p1,p2,p3,p4,p5,p6,align = "h",nrow =3)
plot_grid(p7,p8,p9,p10,align = "h",nrow =2)
plot_grid(p1,p2,p5,p6,align = "h",nrow =2)

COLS <- c("deepskyblue3","tomato3","lightgreen")

library(ggstatsplot)
p7=ggpiestats(rt2,  cluster, Sex, palette = 'Set2',title = 'Gender') +
  scale_fill_manual(values = COLS)
head(rt2)

p8=ggpiestats(rt2,  cluster, ICU, palette = 'Set2',title = 'ICU') + 
  scale_fill_manual(values = COLS) 


head(rt2)

p9=ggpiestats(rt2,  cluster, mechanical_ventilation, palette = 'Set2',title = 'mechanical_ventilation')+ 
  scale_fill_manual(values = COLS)  


library(cowplot)
plot_grid(p7,p8,p9,align = "h",nrow = 3)

