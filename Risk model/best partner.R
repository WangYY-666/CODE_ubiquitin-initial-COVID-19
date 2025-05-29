rm(list = ls())
.libPaths("D:/PackagesR/R_LIBS")
library(dplyr)
load(file = "LASSO_need.Rdata")
set.seed(123) 
# 将尚未在训练集中的行ID分配到验证中
glmnet_input<-housing.df
#glmnet_input<-housing.df
x <- data.matrix(glmnet_input[,3:length(colnames(glmnet_input))])
library(survival)
y <- data.matrix(Surv(glmnet_input[,1],glmnet_input[,2]))
library(glmnet)
library(broom)
library(patchwork) 
library(tidyverse)

roci<-data.frame(a=NA,gene=NA,number=NA,time=NA,AUC=NA)
for (m in 55:100) {
  print(m)
  cv.fit<-cv.glmnet(x, y, family = "cox", alpha = m/100) #建立LASSO回归模型，
  coef1<-coef(cv.fit,s = "lambda.min")
  coef2<-coef(cv.fit,s = "lambda.1se")
  fit <- glmnet(x, y, family = "cox", alpha=m/100)
  fit_predict <- fit
  tidy_df <- broom::tidy(fit)            
  tidy_cvdf <- broom::tidy(cv.fit)            
  tmp <- filter(tidy_df, term != '(Intercept)')            
  if(T){
    (Coefficients <- coef(fit, s = cv.fit$lambda.1se))
    (Active.Index <- which(as.matrix(Coefficients) != 0))
    (Active.Coefficients  <- Coefficients[Active.Index])
    (gene_found <- row.names(Coefficients)[Active.Index])
  }
  if(T){
    coxdata <- cbind(glmnet_input[,c(1,2)],glmnet_input[,gene_found])
    res <- data.frame()
    genes = gene_found
    for (i in 1:length(genes)) {
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
  }
  coxresult <- cbind(res[,],Active.Coefficients)
  colnames(coxresult) <- c("Symbol", 
                           "HR", "95% CI", "p.value", "Cox Coefficient")
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
  }
  
  b<-min(which(survivalROC_data$t == "20")) 
  rocx<-survivalROC_data[c(1,b,nrow(survivalROC_data)),c(1,2)]
  rocx<-cbind(a=c(m/100,NA,NA),gene=c(gene_found[1],gene_found[2],gene_found[3]),number=c(nrow(res),NA,NA),rocx)
  colnames(rocx)<-colnames(roci)
  roci<-rbind(roci,rocx)
}
roci<-roci[-1,]
write.table(roci,file = "robust-patent.txt",sep = "\t",row.names = F,quote = F)  

rt<-data.table::fread("robust-plot.txt")
ggplot(data=rt,aes(x=Group,y=Day10_AUC,colour = Group))+ 
  geom_violin(#color = 'grey',
    alpha = 0.8,#alpha = 0.8 参数控制着小提琴图的透明度。具体来说，这里的 alpha 参数用于指定填充颜色的透明度，数值范围通常在 0 到 1 之间，其中 0 表示完全透明（即不可见），1 表示完全不透明。
    scale = 'width',#小提琴宽度
    #linewidth = 1, #外轮廓粗细
    trim = TRUE)+#trim = TRUE 参数控制着小提琴图的形状。当 trim = TRUE 时，小提琴图会根据数据的分布进行修剪
  geom_boxplot(mapping=aes(x=Group,y=Day10_AUC,colour = Group,fill=Group), #箱线图
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=Day10_AUC,colour = Group), #散点
              alpha = 0.3,size=3)+
  geom_signif(mapping=aes(x=Group,y=Day10_AUC), # 不同组别的显著性
              comparisons = list(c("2 (UBE2V1 ISG15)", "4 (UBE2V1 UBE2L3 UBE2B ISG15)"), # 哪些组进行比较
                                 c("2 (UBE2V1 ISG15)", "3 (UBE2V1 UBE2B ISG15)"),
                                 c("4 (UBE2V1 UBE2L3 UBE2B ISG15)", "3 (UBE2V1 UBE2B ISG15)")),
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = c(0.795,0.775,0.785), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 4, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型,可以更改
  theme_bw()+#设置白色背景
  labs(x="不同基因数",y=" 10天的AUC值") # 添加标题，x轴，y轴内容
ggplot(data=rt,aes(x=Group,y=Day20_AUC,colour = Group))+ 
  geom_violin(#color = 'grey',
    alpha = 0.8,#alpha = 0.8 参数控制着小提琴图的透明度。具体来说，这里的 alpha 参数用于指定填充颜色的透明度，数值范围通常在 0 到 1 之间，其中 0 表示完全透明（即不可见），1 表示完全不透明。
    scale = 'width',#小提琴宽度
    #linewidth = 1, #外轮廓粗细
    trim = TRUE)+#trim = TRUE 参数控制着小提琴图的形状。当 trim = TRUE 时，小提琴图会根据数据的分布进行修剪
  geom_boxplot(mapping=aes(x=Group,y=Day20_AUC,colour = Group,fill=Group), #箱线图
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=Day20_AUC,colour = Group), #散点
              alpha = 0.3,size=3)+
  geom_signif(mapping=aes(x=Group,y=Day20_AUC), # 不同组别的显著性
              comparisons = list(c("2 (UBE2V1 ISG15)", "4 (UBE2V1 UBE2L3 UBE2B ISG15)"), # 哪些组进行比较
                                 c("2 (UBE2V1 ISG15)", "3 (UBE2V1 UBE2B ISG15)"),
                                 c("4 (UBE2V1 UBE2L3 UBE2B ISG15)", "3 (UBE2V1 UBE2B ISG15)")),
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = c(0.76,0.754,0.757), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 4, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型,可以更改
  theme_bw()+#设置白色背景
  labs(x="不同基因数",y=" 20天的AUC值") # 添加标题，x轴，y轴内容
ggplot(data=rt,aes(x=Group,y=Day30_AUC,colour = Group))+ 
  geom_violin(#color = 'grey',
    alpha = 0.8,#alpha = 0.8 参数控制着小提琴图的透明度。具体来说，这里的 alpha 参数用于指定填充颜色的透明度，数值范围通常在 0 到 1 之间，其中 0 表示完全透明（即不可见），1 表示完全不透明。
    scale = 'width',#小提琴宽度
    #linewidth = 1, #外轮廓粗细
    trim = TRUE)+#trim = TRUE 参数控制着小提琴图的形状。当 trim = TRUE 时，小提琴图会根据数据的分布进行修剪
  geom_boxplot(mapping=aes(x=Group,y=Day30_AUC,colour = Group,fill=Group), #箱线图
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=Group,y=Day30_AUC,colour = Group), #散点
              alpha = 0.3,size=3)+
  geom_signif(mapping=aes(x=Group,y=Day30_AUC), # 不同组别的显著性
              comparisons = list(c("2 (UBE2V1 ISG15)", "4 (UBE2V1 UBE2L3 UBE2B ISG15)"), # 哪些组进行比较
                                 c("2 (UBE2V1 ISG15)", "3 (UBE2V1 UBE2B ISG15)"),
                                 c("4 (UBE2V1 UBE2L3 UBE2B ISG15)", "3 (UBE2V1 UBE2B ISG15)")),
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = c(0.86,0.85,0.855), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 4, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型,可以更改
  theme_bw()+#设置白色背景
  labs(x="不同基因数",y=" 30天的AUC值") # 添加标题，x轴，y轴内容
