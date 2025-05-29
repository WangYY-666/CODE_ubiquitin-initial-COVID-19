.libPaths("D:/PackagesR/R_LIBS")
library(ggpubr)
rm(list = ls())
load(file = "risksocre_order.Rdata")
colnames(rt)[5]<-"ID"
rt1<-data.table::fread(file="violin.txt")
colnames(rt1)[1]<-"ID"
rt1<-rt1[,-c(2,6,10)]
library(dplyr)

rt2<-merge(rt1,rt,by="ID")


#刻度朝内（需要改的话，改每个图的后两行，负值改正即可


head(rt2)
# Change color by groups
dp <- ggplot(rt2, aes(x=risk, y=as.numeric(rt2$apacheii), fill=risk)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of APACHEII",x="risk", y = "APACHEII")+ 
  stat_compare_means() # Add global p-value
p1=dp + scale_fill_manual(values=c("high"="#D95F02","low"="#1B9E77"))+
  scale_y_continuous(limits = c(0,55))+ theme_minimal()+ 
  #theme(axis.line=element_line(colour="black"))+
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 0,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")+
  theme(axis.ticks.length.x=unit(-0.25,'cm'))+
  theme(axis.ticks.length.y=unit(-0.25,'cm'))
p1

head(rt2)
# Change color by groups
dp <- ggplot(rt2, aes(x=risk, y=as.numeric(rt2$sofa), fill=risk)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of SOFA",x="risk", y = "SOFA")+ 
  stat_compare_means() # Add global p-value
p2=dp + scale_fill_manual(values=c("high"="#D95F02","low"="#1B9E77"))+
  scale_y_continuous(limits = c(0,25))+ theme_minimal()+ 
  #theme(axis.line=element_line(colour="black"))+
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 0,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")+
  theme(axis.ticks.length.x=unit(-0.25,'cm'))+
  theme(axis.ticks.length.y=unit(-0.25,'cm'))

p2
head(rt2)
# Change color by groups
dp <- ggplot(rt2, aes(x=risk, y=as.numeric(rt2$HFDP45), fill=risk)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of HFD45",x="risk", y = "HFD45")+ 
  stat_compare_means() # Add global p-value
p3=dp + scale_fill_manual(values=c("high"="#D95F02","low"="#1B9E77"))+
  scale_y_continuous(limits = c(0,60))+ theme_minimal()+ 
  #theme(axis.line=element_line(colour="black"))+
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 0,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")+
  theme(axis.ticks.length.x=unit(-0.25,'cm'))+
  theme(axis.ticks.length.y=unit(-0.25,'cm'))

p3
head(rt2)
# Change color by groups
dp <- ggplot(rt2, aes(x=risk, y=as.numeric(rt2$ventilator_free_days), fill=risk)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of ventilator_free_days",x="risk", y = "ventilator_free_days")+ 
  stat_compare_means() # Add global p-value
p4=dp + scale_fill_manual(values=c("high"="#D95F02","low"="#1B9E77"))+
  scale_y_continuous(limits = c(0,45))+ theme_minimal()+ 
  #theme(axis.line=element_line(colour="black"))+
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 0,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")+
  theme(axis.ticks.length.x=unit(-0.25,'cm'))+
  theme(axis.ticks.length.y=unit(-0.25,'cm'))


p4


head(rt2)
# Change color by groups
dp <- ggplot(rt2, aes(x=risk, y=as.numeric(rt2$`ferritin(ng/ml)`), fill=risk)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of Ferritin",x="risk", y = "Ferritin(ng/ml)")+ 
  stat_compare_means() # Add global p-value
p5=dp + scale_fill_manual(values=c("high"="#D95F02","low"="#1B9E77"))+ theme_minimal()+ 
  #theme(axis.line=element_line(colour="black"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits = c(0,6000))+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 0,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")+
  theme(axis.ticks.length.x=unit(-0.25,'cm'))+
  theme(axis.ticks.length.y=unit(-0.25,'cm'))

p5

head(rt2)
# Change color by groups
dp <- ggplot(rt2, aes(x=risk, y=as.numeric(rt2$`crp (mg/l)`), fill=risk)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of CRP",x="risk", y = "CRP (mg/l)")+ 
  stat_compare_means() # Add global p-value
p6=dp + scale_fill_manual(values=c("high"="#D95F02","low"="#1B9E77"))+ theme_minimal()+ 
  #theme(axis.line=element_line(colour="black"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits = c(0,550))+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 0,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")+
  theme(axis.ticks.length.x=unit(-0.25,'cm'))+
  theme(axis.ticks.length.y=unit(-0.25,'cm'))
p6

head(rt2)
# Change color by groups
dp <- ggplot(rt2, aes(x=risk, y=as.numeric(rt2$`ddimer (mg/l_feu)`), fill=risk)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.01, fill="white")+
  labs(title="Plot of DDIMER",x="risk", y = "DDIMER (mg/l_feu)")+ 
  stat_compare_means() # Add global p-value
p7=dp + scale_fill_manual(values=c("high"="#D95F02","low"="#1B9E77"))+
  scale_y_continuous(limits = c(0,120))+ theme_minimal()+ 
  #theme(axis.line=element_line(colour="black"))+
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 0,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")+
  theme(axis.ticks.length.x=unit(-0.25,'cm'))+
  theme(axis.ticks.length.y=unit(-0.25,'cm'))

p7
head(rt2)
# Change color by groups
dp <- ggplot(rt2, aes(x=risk, y=as.numeric(rt2$`procalcitonin (ng/ml)`), fill=risk)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of Procalcitonin",x="risk", y = "Procalcitonin (ng/ml)")+ 
  stat_compare_means() # Add global p-value
p8=dp + scale_fill_manual(values=c("high"="#D95F02","low"="#1B9E77"))+ theme_minimal()+ 
  #theme(axis.line=element_line(colour="black"))+
  theme(axis.title.x = element_blank())+theme_classic()+
  theme(axis.title.y = element_text(size = 15,vjust = 0.5,hjust = 0.5))+
  theme(axis.text.x = element_text(size = 15,angle = 0,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(legend.position = "NA")+
  theme(axis.ticks.length.x=unit(-0.25,'cm'))+
  theme(axis.ticks.length.y=unit(-0.25,'cm'))



p8

library(cowplot)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,align = "h",nrow =2)




library(ggstatsplot)
rt2[is.na(rt2)]=0

head(rt2)
p9=ggbarstats(rt2, risk, ICU, palette = 'Set2')

head(rt2)
p10=ggbarstats(rt2, risk,mechanical_ventilation, palette = 'Set2')

load(file = "clust_MR.Rdata")
rt3<-as.data.frame(clust_mr)
colnames(rt3)<-c("ID","cluster")
rt4<-merge(rt3,rt2,by="ID")
table(rt4$cluster)
p11=ggbarstats(rt4, risk,cluster, palette = 'Set2')

library(cowplot)
plot_grid(p9,p10,p11,align = "h",nrow =1)
#plot_grid(p9,p10,align = "h",nrow =1)
