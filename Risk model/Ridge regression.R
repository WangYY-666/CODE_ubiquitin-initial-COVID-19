rm(list = ls())
.libPaths("D:/PackagesR/R_LIBS")

library(dplyr)
#load(file = "LASSO_need.Rdata")
load(file = "GSE157103_survivaltest_normalized.Rdata")
gene<-c("UBE2V1",  "UBE2L3" ,"UBE2B" ,"ISG15"  )
#gene<-c("UBE2V1","UBE2B" )
rt_cli<-expr
rt_cli1<-rt_cli[,gene]

rt_cli2<-cbind(rt_cli[,c(2,3)],rt_cli1)
housing.df<-as.data.frame(lapply(rt_cli2,as.numeric))
rownames(housing.df)<-rt_cli$id

#housing.df<-expr[,2:19475]
#rownames(housing.df)<-expr[,1]
save(housing.df,file = "LASSO_need.Rdata")
#改
housing.df<-housing.df%>%
  select(2,1,3:6)

set.seed(123) ## to get the same sequence of numbers
# randomly sample 60% of the row IDs for training; the remaining 40% serve as validation
train.rows <- sample(rownames(housing.df), dim(housing.df)[1]*0.55)
housing.train <- housing.df[train.rows, ]
# 将尚未在训练集中的行ID分配到验证中
valid.rows <- setdiff(rownames(housing.df), train.rows)
housing.valid <- housing.df[valid.rows, ]
save(train.rows,file = "train_rows.Rdata")
glmnet_input<-housing.train
#glmnet_input<-housing.df
x <- data.matrix(glmnet_input[,3:length(colnames(glmnet_input))])

library(survival)
y <- data.matrix(Surv(glmnet_input[,1],glmnet_input[,2]))
## lasso
library(glmnet)
#cv.fit <- cv.glmnet(x, y, family="cox") 
cv.fit<-cv.glmnet(x, y, family = "cox", alpha = 0.33) #建立LASSO回归模型，

cv.fit$lambda.min
cv.fit$lambda.1se

#x11()
plot(cv.fit)

if(T){
  library(ggplot2)
  library(broom)
  tidied_cv = tidy(cv.fit)
  ggplot(tidied_cv, aes(lambda, estimate)) +
    geom_point(color = "red",size=2) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), alpha = .3) +
    scale_x_log10()+
    geom_vline(xintercept = cv.fit$lambda.1se, lty=2,color = "red",size=1) +
    geom_vline(xintercept = cv.fit$lambda.min, lty=2) +
    ylab("Partial Likelihood Deviance")+
    annotate("text",x=cv.fit$lambda.1se/5,size=5,
             y =max(tidied_cv$estimate),
             label =paste0("lambda.1se = ",
                           round(cv.fit$lambda.1se,3),
                           "\n",
                           "number = ",tidied_cv$nzero[tidied_cv$lambda==cv.fit$lambda.1se]))+
    ggtitle("")+
    theme_bw()
}

ggsave("图1-train-Cross validation curve.pdf",width =2,height = 2)
coef1<-coef(cv.fit,s = "lambda.min")
coef2<-coef(cv.fit,s = "lambda.1se")
coef1
coef2
fit <- glmnet(x, y, family = "cox", alpha=0.33)

fit_predict <- fit

library(broom)
library(patchwork) 
library(tidyverse)
tidy_df <- broom::tidy(fit)            
tidy_cvdf <- broom::tidy(cv.fit)            
tmp <- filter(tidy_df, term != '(Intercept)')            
ggplot(tmp, aes(log(lambda),estimate,color=term,group=term))+     geom_line(size=1)+            
  labs(x="Log Lambda",y="Coefficients")+            
  coord_cartesian(xlim = c(-9,-1), expand = T) +            
  guides(color = guide_legend(title = ''),direction = 'cox') +            
  #scale_x_continuous(limits = c(-9,-1)) +            
  theme_bw()

plot(fit, label = TRUE)
if(T){
  (Coefficients <- coef(fit, s = cv.fit$lambda.1se))
  (Active.Index <- which(as.matrix(Coefficients) != 0))
  (Active.Coefficients  <- Coefficients[Active.Index])
  
  (gene_found <- row.names(Coefficients)[Active.Index])
  save(gene_found,file = "gene_found.Rdata")
}

if(T){
  tidy_covariates = tidy(cv.fit$glmnet.fit)
  index <- tidy_covariates$term %in% gene_found
  library(dplyr)
  data_point <- tidy_covariates %>% 
    filter(index) %>% 
    arrange(lambda) %>% 
    distinct(term,.keep_all = T)
  library(ggrepel)
  library(ggplot2)
  ggplot(tidy_covariates)+
    geom_line(data = subset(tidy_covariates,index), aes(x=lambda, y=estimate, color=as.factor(term)),size=1) +
    geom_line(data = subset(tidy_covariates,!index), aes(x=lambda, y=estimate, color=as.factor(term)),size=0.8,alpha=0.1) +
    guides(color=FALSE) +
    geom_vline(xintercept = cv.fit$lambda.1se, lty=2,color='red') +
    scale_x_log10()+
    ylab("Coefficients")+
    geom_point(data =data_point,aes(lambda,estimate,color=as.factor(term)),size=3)+
    geom_text_repel(data =data_point,aes(lambda,estimate,label=term),size=3)+
    theme_classic()
}
ggsave("图2-train-Coefficients.pdf",width =2,height = 2)
if(T){
  coxdata <- cbind(glmnet_input[,c(1,2)],glmnet_input[,gene_found])
  res <- data.frame()
  genes = gene_found
  for (i in 1:length(genes)) {
    print(i)
    surv =as.formula(paste('Surv(ventilator_free_days,mechanical_ventilation)~', genes[i]))
    cox = coxph(surv, data = coxdata)
    cox = summary(cox)
    p.value=signif(cox$wald["pvalue"], digits=2)
    HR =signif(cox$coef[2], digits=2);#exp(beta)
    HR.confint.lower = signif(cox$conf.int[,"lower .95"], 2)
    HR.confint.upper = signif(cox$conf.int[,"upper .95"],2)
    CI <- paste0("(", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res[i,1] = genes[i]
    res[i,2] = HR
    res[i,3] = CI
    res[i,4] = p.value
  }
  names(res) <- c("ID","HR","95% CI","p.value")
  save(res,file = "res_HR.Rdata")
}


coxresult <- cbind(res[,],Active.Coefficients)

colnames(coxresult) <- c("Symbol", 
                         "HR", "95% CI", "p.value", "Cox Coefficient")
write.table(coxresult,file = "train-coefficience.txt",sep = "\t",row.names = F,quote = F)
coxresult$HR <- as.numeric(as.character(coxresult$HR))
coxresult <-  dplyr::arrange(coxresult,desc(HR))


x <- data.matrix(glmnet_input[,3:length(colnames(glmnet_input))])

rt <- glmnet_input  
rt <- rt[,1:2]
order =rownames(rt)
library(glmnet)
rt$riskScore <- predict(fit_predict, x, s=cv.fit$lambda.1se, type="link")
# best AUC cutoff
rt <- as.data.frame(apply(rt,2,as.numeric))
rownames(rt) <- rownames(glmnet_input)

# rt1<-cbind(id=rownames(rt),rt)
# rt2<-cbind(id=rownames(glmnet_input),glmnet_input[,c(3:6)])
#rt3<-merge(rt1,rt2,by="id")
#save(rt3,file = "cox_riskscore.Rdata")


source("survivalROC_NEW.R")
cutoffROC <- function(t,rt){
  rt=rt
  ROC_rt <- survivalROC_NEW(Stime=rt[,1], status=rt[,2], marker = rt[,3], 
                            predict.time =t, method="KM")
  
  youdenindex=max(ROC_rt$sensitivity+ROC_rt$specificity-1)
  
  index = which(ROC_rt$sensitivity+ROC_rt$specificity-1==youdenindex)+1
  
  result=c(ROC_rt$cut.values[index],ROC_rt$sensitivity[index-1],ROC_rt$specificity[index-1])
}

t <- 25#平均值
(res.cut <- cutoffROC(t,rt))
cutoff <- res.cut[1]



library(dplyr)
if(T){
  rt$risk=as.vector(ifelse(rt$riskScore>cutoff,"high","low"))
  rt$order <- order
  rt$riskScore = as.numeric(rt$riskScore)
  rt = rt %>% arrange(riskScore)
  rt$num =seq(1,length(rownames(rt)))
}
rt1<-rt
####
#save(rt1,file = "18need.Rdata") 

if(T){
  cutnum <- sum(rt$risk=="low")
  library(ggplot2)
  (plot.point=ggplot(rt,aes(x=num,y=riskScore))+
      geom_point(aes(col=risk),size=3)+
      geom_segment(aes(x = cutnum, y = min(riskScore), xend = cutnum, yend = cutoff),linetype="dashed")+
      geom_segment(aes(x = 0, y = cutoff, xend = cutnum,yend = cutoff),linetype="dashed")+
      geom_text(aes(x=cutnum/2,y=cutoff+0.25,label=paste0("Cutoff: ",round(cutoff,2))),col ="black",size = 8,alpha=0.8)+
      theme_classic()+
      theme(axis.title.x=element_blank(),#legend.text = element_text(size = 15),
            axis.text.x=element_text(size = 20),
            #axis.ticks.x=element_blank(),
            axis.title.y =  element_text(size = 20),
            axis.text.y = element_text(size = 15))+
      scale_color_manual(values=c("#FC4E07","#00AFBB"))+
      scale_x_continuous(limits = c(0,NA),expand = c(0,0))
  )
}


if(T){
  event = factor(rt$mechanical_ventilation)
  (plot.sur=ggplot(rt,aes(x=num,y=ventilator_free_days))+
      geom_point(aes(col=event),size=3)+
      geom_vline(aes(xintercept = cutnum),linetype="dashed")+
      theme_classic()+
      theme(axis.title.x=element_text(size = 15),
            axis.text.x=element_text(size = 20),
            axis.ticks.x=element_blank(),
            axis.title.y =  element_text(size = 20),
            axis.text.y = element_text(size = 15))+      
      scale_color_manual(values=c("#00AFBB","#FC4E07"))+
      scale_x_continuous(limits = c(0,NA),expand = c(0,0)))
}


if(T){
  tmp=t(glmnet_input[rt$order,rev(coxresult$Symbol)])
  library(data.table)
  library(plyr)
  library(scales)
  tmp.m <- reshape2::melt(tmp)
  colnames(tmp.m)=c("Name", "variable", "value")
  tmp.m <- ddply(tmp.m, .(variable), transform,rescale = rescale(value))
  (plot.h <- ggplot(tmp.m, aes(variable, Name)) + 
      geom_tile(aes(fill = rescale), colour = "white") + 
      scale_fill_gradient(low = "#00AFBB",high = "#FC4E07")+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y = element_text(size = 15))
  )
}

library(cowplot)
plot_grid(plot.point, plot.sur, plot.h,
          align = 'v',ncol = 1,axis="b")
ggsave("图3-train-riskscore.pdf",width =8,height = 10)
#ggsave("图3-riskscore.pdf",width =8,height = 10)
##########################################################################################



if(T){
  library("survival")
  library("survminer")
  my.surv <- Surv(rt$ventilator_free_days, rt$mechanical_ventilation)
  group <- rt$risk
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  
  
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  x = summary(coxph(Surv(ventilator_free_days, mechanical_ventilation)~riskScore, data = rt))
  HR = signif(x$coef[2], digits=2)
  up95 = signif(x$conf.int[,"upper .95"],2)
  low95 = signif(x$conf.int[,"lower .95"], 2)
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  rt <- rt[order(rt[,"riskScore"],decreasing = T),]
  plot.km<-ggsurvplot(fit, data = survival_dat ,
             
             conf.int = F, 
             
             censor = F,
             palette = c("#D95F02","#1B9E77"), 
             
             font.legend = 11,
             
             legend.labs=c(paste0(">",round(cutoff,2),"(",sum(rt$risk=="high"),")"),
                           paste0("<",round(cutoff,2),"(",sum(rt$risk=="low"),")")),
             ylab=c('Mechanical ventilation probability '),
             
             pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                        paste("p = ",round(p.val,3), sep = "")),
                          HR, CI, sep = "\n"))
}
ggsave("图4-train-生存分析.pdf",width = 5.5,height = 5.5)
#ggsave("图4-生存分析.pdf",width = 5.5,height = 5.5)

##ROC
if(T){
  library(tidyr)
  library(purrr)
  library(survivalROC)
  ## Define a helper function to evaluate at various t
  survivalROC_helper <- function(t) {
    survivalROC(Stime=rt$ventilator_free_days, status=rt$mechanical_ventilation, marker = rt$riskScore, 
                predict.time =t, method="KM")
  }
  ## Evaluate 3,5,10
  survivalROC_data <- data_frame(t = c(10,20,30)) %>%
    mutate(survivalROC = map(t, survivalROC_helper),
           ## Extract scalar AUC
           auc = purrr::map_dbl(survivalROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_survivalROC = map(survivalROC, function(obj) {
             dplyr::as_data_frame(obj[c("cut.values","TP","FP")])
           })) %>%
    dplyr::select(-survivalROC) %>%
    unnest() %>%
    arrange(t, FP, TP)
  ## Plot
  survivalROC_data1 <- survivalROC_data %>% 
    mutate(auc =sprintf("%.3f",auc))%>% 
    unite(year, t,auc,sep = " days AUC: ")
  
  AUC =factor(survivalROC_data1$year)
  plot.roc<-survivalROC_data1 %>%
    ggplot(mapping = aes(x = FP, y = TP)) +
    geom_path(aes(color= AUC),size=2)+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    theme_classic() +
    theme(legend.position = c(0.7,0.2))+theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank(),
                                              axis.title.y =  element_text(size = 15),
                                              axis.text.y = element_text(size = 15),
                                              legend.text = element_text(size = 12))
}

ggsave("图5-train-ROC.pdf",width = 5.5,height = 5.0)
#ggsave("图5-ROC.pdf",width = 5.5,height = 5.0)
plot_grid(plot.km$plot, plot.roc,
          align = 'v',ncol = 1,axis="b")
ggsave("图6-train-sum.pdf",width =8,height = 10)











#***validation
glmnet_input=housing.valid
x <- data.matrix(glmnet_input[,3:length(colnames(glmnet_input))])

rt <- glmnet_input  
rt <- rt[,1:2]
order =rownames(rt)

library(glmnet)
rt$riskScore <- predict(fit_predict, x, s=cv.fit$lambda.1se, type="link")

rt <- as.data.frame(apply(rt,2,as.numeric))
rownames(rt) <- rownames(glmnet_input)

if(T){
  rt$risk=as.vector(ifelse(rt$riskScore>cutoff,"high","low"))
  rt$order <- order
  rt$riskScore = as.numeric(rt$riskScore)
  rt = rt %>% arrange(riskScore)
  rt$num =seq(1,length(rownames(rt)))
}
rt2<-rt
####

source("survivalROC_NEW.R")



# 准备数据
library(dplyr)
if(T){
  rt$risk=as.vector(ifelse(rt$riskScore>cutoff,"high","low"))
  rt$order <- order
  rt$riskScore = as.numeric(rt$riskScore)
  rt = rt %>% arrange(riskScore)
  rt$num =seq(1,length(rownames(rt)))
}


if(T){
  cutnum <- sum(rt$risk=="low")
  library(ggplot2)
  (plot.point=ggplot(rt,aes(x=num,y=riskScore))+
      geom_point(aes(col=risk),size=3)+
      geom_segment(aes(x = cutnum, y = min(riskScore), xend = cutnum, yend = cutoff),linetype="dashed")+
      geom_segment(aes(x = 0, y = cutoff, xend = cutnum,yend = cutoff),linetype="dashed")+
      geom_text(aes(x=cutnum/2,y=cutoff+0.25,label=paste0("Cutoff: ",round(cutoff,2))),col ="black",size =8,alpha=0.8)+
      theme_classic()+
      theme(axis.title.x=element_blank(),#legend.text = element_text(size = 15),
            axis.text.x=element_text(size = 20),
            #axis.ticks.x=element_blank(),
            axis.title.y =  element_text(size = 20),
            axis.text.y = element_text(size = 15))+
      scale_color_manual(values=c("#FC4E07","#00AFBB"))+
      scale_x_continuous(limits = c(0,NA),expand = c(0,0)))
}


if(T){
  event = factor(rt$mechanical_ventilation)
  (plot.sur=ggplot(rt,aes(x=num,y=ventilator_free_days))+
      geom_point(aes(col=event),size=3)+
      geom_vline(aes(xintercept = cutnum),linetype="dashed")+
      theme_classic()+
      theme(axis.title.x=element_text(size = 15),
            axis.text.x=element_text(size = 20),
            axis.ticks.x=element_blank(),
            axis.title.y =  element_text(size = 20),
            axis.text.y = element_text(size = 15))+
      scale_color_manual(values=c("#00AFBB","#FC4E07"))+
      scale_x_continuous(limits = c(0,NA),expand = c(0,0)))
}


if(T){
  tmp=t(glmnet_input[rt$order,rev(coxresult$Symbol)])
  library(data.table)
  library(plyr)
  library(scales)
  tmp.m <- reshape2::melt(tmp)
  colnames(tmp.m)=c("Name", "variable", "value")
  tmp.m <- ddply(tmp.m, .(variable), transform,rescale = rescale(value))
  (plot.h <- ggplot(tmp.m, aes(variable, Name)) + 
      geom_tile(aes(fill = rescale), colour = "white") + 
      scale_fill_gradient(low = "#00AFBB",high = "#FC4E07")+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y = element_text(size = 15)))
}

library(cowplot)
plot_grid(plot.point, plot.sur, plot.h,
          
          align = 'v',ncol = 1,axis="b")
ggsave("图7-validate-riskscores.pdf",width =8,height = 10)
##########################################################################################



if(T){
  library("survival")
  library("survminer")
  my.surv <- Surv(rt$ventilator_free_days, rt$mechanical_ventilation)
  group <- rt$risk
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  
  
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  x = summary(coxph(Surv(ventilator_free_days, mechanical_ventilation)~riskScore, data = rt))
  HR = signif(x$coef[2], digits=2)
  up95 = signif(x$conf.int[,"upper .95"],2)
  low95 = signif(x$conf.int[,"lower .95"], 2)
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  
  rt <- rt[order(rt[,"riskScore"],decreasing = T),]
  plot.km<-ggsurvplot(fit, data = survival_dat ,
             
             conf.int = F, 
             
             censor = F, 
             palette = c("#D95F02","#1B9E77"), 
             
             font.legend = 11,
             
             legend.labs=c(paste0(">",round(cutoff,2),"(",sum(rt$risk=="high"),")"),
                           paste0("<",round(cutoff,2),"(",sum(rt$risk=="low"),")")),
             ylab=c('mechanical ventilation probability '),
             
             pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                        paste("p = ",round(p.val,3), sep = "")),
                          HR, CI, sep = "\n"))
}
ggsave("图8-validate-生存分析.pdf",width = 5.5,height = 5.5)

##ROC
if(T){
  library(tidyr)
  library(purrr)
  library(survivalROC)
  ## Define a helper function to evaluate at various t
  survivalROC_helper <- function(t) {
    survivalROC(Stime=rt$ventilator_free_days, status=rt$mechanical_ventilation, marker = rt$riskScore, 
                predict.time =t, method="KM")
  }
  ## Evaluate 3,5,10
  survivalROC_data <- data_frame(t = c(10,20,30)) %>%
    mutate(survivalROC = map(t, survivalROC_helper),
           ## Extract scalar AUC
           auc = purrr::map_dbl(survivalROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_survivalROC = map(survivalROC, function(obj) {
             dplyr::as_data_frame(obj[c("cut.values","TP","FP")])
           })) %>%
    dplyr::select(-survivalROC) %>%
    unnest() %>%
    arrange(t, FP, TP)
  ## Plot
  survivalROC_data1 <- survivalROC_data %>% 
    mutate(auc =sprintf("%.3f",auc))%>% 
    unite(year, t,auc,sep = " days AUC: ")
  
  AUC =factor(survivalROC_data1$year)
  plot.roc<-survivalROC_data1 %>%
    ggplot(mapping = aes(x = FP, y = TP)) +
    geom_path(aes(color= AUC),size=2)+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    theme_classic() +
    theme(legend.position = c(0.7,0.2))+theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank(),
                                              axis.title.y =  element_text(size = 15),
                                              axis.text.y = element_text(size = 15),
                                              legend.text = element_text(size = 12))
}

ggsave("图9-validate-ROC.pdf",width = 5.5,height = 5.0)
plot_grid(plot.km$plot, plot.roc,
          align = 'v',ncol = 1,axis="b")
ggsave("图10-validate-sum.pdf",width =8,height = 10)
#save(res,file = "symbolgene.Rdata")

###预测模型评分比分成高低两组
glmnet_input<-housing.df
x <- data.matrix(glmnet_input[,3:length(colnames(glmnet_input))])

rt <- glmnet_input  
rt <- rt[,1:2]
order =rownames(rt)

library(glmnet)
rt$riskScore <- predict(fit_predict, x, s=cv.fit$lambda.1se, type="link")

rt <- as.data.frame(apply(rt,2,as.numeric))
rownames(rt) <- rownames(glmnet_input)
if(T){
  rt$risk=as.vector(ifelse(rt$riskScore>cutoff,"high","low"))
  rt$order <- order
  rt$riskScore = as.numeric(rt$riskScore)
  rt = rt %>% arrange(riskScore)
  rt$num =seq(1,length(rownames(rt)))
}
save(rt,file = "risksocre_order.Rdata")










