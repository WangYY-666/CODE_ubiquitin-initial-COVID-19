#https://mp.weixin.qq.com/s/y7ejxowWt9yBwZtE0Z7Tzw
.libPaths("D:/PackagesR/R_LIBS")
rm(list = ls())

load(file = "GSE157103_survivaltest_normalized.Rdata")
load(file = "clust_MR.Rdata")
#cli<-data.table::fread("violin.txt")
#colnames(cli)[1]<-"id"
#x<-merge(clust_mr,cli[,c(1,7),],by="id")
#y<-merge(x,expr[,c(1:3)],by="id")
load(file = "GSE157103_Clinic_tidy.Rdata")
rtt$mechanical_ventilation[rtt$mechanical_ventilation == "yes"] <- 1 
rtt$mechanical_ventilation[rtt$mechanical_ventilation == "no"] <- 0 
rtt$ventilator_free_days<-as.numeric(rtt$ventilator_free_days)
rtt$mechanical_ventilation<-as.numeric(rtt$mechanical_ventilation)
rt<-merge(clust_mr,rtt[,c(1,7,8)],by="id")
rt<-merge(clust_mr,expr[,c(1:3)],by="id")

test<-data.table::fread("GSE157103_clinic_new.txt")
test<-test[,c(1,2,8)]
for (i in 1:nrow(test)) {
  if(test$HFDP45[i]==0 || test$HFDP45[i]>45)
    test$status[i]<-0
  else
    test$status[i]<-1
}
rt<-test
test<-data.table::fread("GSE157103_clinic_new.txt")
test<-test[,c(1,2,7,11)]
for (i in 1:nrow(test)) {
  if(test$mechanical_ventilation[i]=="Yes")
    test$mechanical_ventilation[i]<- 1
  if(test$mechanical_ventilation[i]=="No")
    test$mechanical_ventilation[i]<- 0
}
test$mechanical_ventilation<-as.numeric(test$mechanical_ventilation)
rt<-test
for (i in 1:nrow(rt)) {
  if(rt[i,2]=="1")
    rt[i,2]<-"cluster1"
  if(rt[i,2]=="2")
    rt[i,2]<-"cluster2"
  if(rt[i,2]=="3")
    rt[i,2]<-"cluster3"
}
colnames(rt)[2]<-"cluster"
library("survival")
library("survminer")

#fit <- survfit(Surv(HFDP45, status) ~ cluster, data = rt)
fit <- survfit(Surv(ventilator_free_days, mechanical_ventilation) ~ cluster, data = rt)
print(fit)
# Summary of survival curves
summary(fit)
# Access to the sort summary table
summary(fit)$table
d <- data.frame(time = fit$time,
                n.risk = fit$n.risk,
                n.event = fit$n.event,
                n.censor = fit$n.censor,
                surv = fit$surv,
                upper = fit$upper,
                lower = fit$lower
)
# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           xlab = "Time(days)",
           ylab = "Ventilator free probability",
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), 
           palette = c("lightgreen","red2", "#2E9FDF"))
