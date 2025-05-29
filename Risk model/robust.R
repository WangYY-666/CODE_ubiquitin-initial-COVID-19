rm(list = ls())
.libPaths("D:/PackagesR/R_LIBS")
library(ggplot2)
library(dplyr)
load(file = "LASSO_need.Rdata")
set.seed(123) 
roci<-data.frame(n=NA,time=NA,AUC=NA)
survivalROCx<-data.frame(n=NA,year=NA,cut.values=NA,TP=NA,FP=NA)
for (i in 1:500) {
  #load(file = "LASSO结果/seed123-alpha0.33/train55-valid45/train_rows.Rdata")
  train.rows <- sample(rownames(housing.df), dim(housing.df)[1]*0.55)
  housing.train <- housing.df[train.rows, ]
  #valid.rows <- setdiff(rownames(housing.df), train.rows)
  #housing.valid <- housing.df[valid.rows, ]
  glmnet_input=housing.train
  rt <- glmnet_input  
  rt <- rt[,1:2]
  order =rownames(rt)
  
  library(glmnet)
  rt$riskScore <- (-0.513275372123213)*glmnet_input$UBE2V1+0.303704288919051*glmnet_input$UBE2B+(-0.209083600450345)*glmnet_input$ISG15
  rt <- as.data.frame(apply(rt,2,as.numeric))
  rownames(rt) <- rownames(glmnet_input)
  cutoff<- -3.33
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
    
    AUC =factor(survivalROC_data1$year)
    survivalROC_data1 %>%
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
  survivalROC_data1<-cbind(n=rep(i,nrow(survivalROC_data1)),survivalROC_data1)
  survivalROCx<-rbind(survivalROCx,survivalROC_data1)
  rocx<-survivalROC_data[c(1,46,91),c(1,2)]
  rocx<-cbind(n=c(i,i,i),rocx)
  colnames(rocx)<-colnames(roci)
  roci<-rbind(roci,rocx)
  
}
roc<-roci[-1,]
survivalROCx<-survivalROCx[-1,]

save(roc,file = "LASSO结果/seed123-alpha0.33/train55-valid45/ROC.Rdata")
save(survivalROCx,file = "LASSO结果/seed123-alpha0.33/train55-valid45/survivalROCx.Rdata")
write.table(roc,file ="LASSO结果/seed123-alpha0.33/train55-valid45/ROC1.txt",sep = "\t",row.names =  FALSE)


roci<-data.frame(proportion=NA,n=NA,time=NA,AUC=NA)
for (j in 45:100) {
  for (i in 1:50) {
    #load(file = "LASSO结果/seed123-alpha0.33/train55-valid45/train_rows.Rdata")
    train.rows <- sample(rownames(housing.df), dim(housing.df)[1]*(j/100))
    housing.train <- housing.df[train.rows, ]
    #valid.rows <- setdiff(rownames(housing.df), train.rows)
    #housing.valid <- housing.df[valid.rows, ]
    glmnet_input=housing.train
    rt <- glmnet_input  
    rt <- rt[,1:2]
    order =rownames(rt)
    
    library(glmnet)
    rt$riskScore <- (-0.513275372123213)*glmnet_input$UBE2V1+0.303704288919051*glmnet_input$UBE2B+(-0.209083600450345)*glmnet_input$ISG15
    rt <- as.data.frame(apply(rt,2,as.numeric))
    rownames(rt) <- rownames(glmnet_input)
    cutoff<- -3.33
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
      
      AUC =factor(survivalROC_data1$year)
      survivalROC_data1 %>%
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
    b<-min(which(survivalROC_data$t == "20")) 
    rocx<-survivalROC_data[c(1,b,nrow(survivalROC_data)),c(1,2)]
    rocx<-cbind(proportion=j/100,n=c(i,i,i),rocx)
    colnames(rocx)<-colnames(roci)
    roci<-rbind(roci,rocx)
    
  }
  
}
roc2<-roci[-1,]
save(roc2,file = "LASSO结果/seed123-alpha0.33/train55-valid45/ROC2.Rdata")
write.table(roc2,file ="LASSO结果/seed123-alpha0.33/train55-valid45/ROC2.txt",sep = "\t",row.names =  FALSE)


rm(list = ls())
load(file = "LASSO结果/seed123-alpha0.33/train55-valid45/ROC.Rdata")
load(file = "LASSO结果/seed123-alpha0.33/train55-valid45/survivalROCx.Rdata")
library(tidyr)  
library(dplyr) 
survivalrt<-separate(survivalROCx, col = "year", into = c("time", "AUC"), sep = ":") 
survivalrt10<-survivalrt[survivalrt$time=="10 days AUC",]
survivalrt20<-survivalrt[survivalrt$time=="20 days AUC",]
survivalrt30<-survivalrt[survivalrt$time=="30 days AUC",]
colours_for_ROC_curves<-rainbow(n=500)
ggplot(mapping = aes(x = FP, y = TP),data=survivalrt10) +
  geom_path(aes(color= n),size=2)+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_classic() +
  theme(legend.position = c(0.7,0.2))+theme(axis.title.x=element_blank(),
                                            axis.text.x=element_blank(),
                                            axis.ticks.x=element_blank(),
                                            axis.title.y =  element_text(size = 15),
                                            axis.text.y = element_text(size = 15),
                                            legend.text = element_text(size = 12))
roc1_ <- roc %>%  
  pivot_wider(  
    id_cols = n,         # 保持id列不变  
    names_from = time,    # 使用name列的值作为新列名  
    values_from = AUC   # 使用value列的值填充新列  
  )  
roc1_<- unlist(sapply(roc1_, as.numeric))
colnames(roc1_)<-c( 'Order',"Day10_AUC","Day20_AUC","Day30_AUC")

roc1_<-as.data.frame(roc1_)

data<-roc[,-1]
colnames(data)<-c("group","value")
data$value<-as.numeric(data$value)
data$group<-as.character(data$group)
ggplot(data=data,aes(x=group,y=value,colour = group))+ 
  geom_violin(#color = 'grey',
    alpha = 0.8,#alpha = 0.8 参数控制着小提琴图的透明度。具体来说，这里的 alpha 参数用于指定填充颜色的透明度，数值范围通常在 0 到 1 之间，其中 0 表示完全透明（即不可见），1 表示完全不透明。
    scale = 'width',#小提琴宽度
    #linewidth = 1, #外轮廓粗细
    trim = TRUE)+#trim = TRUE 参数控制着小提琴图的形状。当 trim = TRUE 时，小提琴图会根据数据的分布进行修剪
  geom_boxplot(mapping=aes(x=group,y=value,colour=group,fill=group), #箱线图
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=group,y=value,colour = group), #散点
              alpha = 0.3,size=3)+
  geom_signif(mapping=aes(x=group,y=value), # 不同组别的显著性
              comparisons = list(c("10", "20"), # 哪些组进行比较
                                 c("10", "30"),
                                 c("20", "30")),
              map_signif_level=T, # T显示显著性，F显示p value
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # 修改显著性线两端的长短
              y_position = c(1.2,1.1,1), # 设置显著性线的位置高度
              size=1, # 修改线的粗细
              textsize = 4, # 修改显著性标记的大小
              test = "t.test")+ # 检验的类型,可以更改
  theme_bw()+#设置白色背景
  labs(x="Days",y="AUC value") # 添加标题，x轴，y轴内容
library(devtools)
install_github("zzawadz/customLayout")
library(customLayout)
lay<- lay_new(mat = matrix(1:3, ncol = 1, nrow = 3),widths = c(1),heights = c(1, 1,1))
lay_show(lay)
lay_set(lay)# 设定绘图对象布局
A<-hist(roc1_$Day10_AUC, prob = T) #直方图
lines(density(na.omit(roc1_$Day10_AUC)), col = 'lightcoral') 
B<-hist(roc1_$Day20_AUC, prob = T) #直方图
lines(density(na.omit(roc1_$Day20_AUC)), col = 'springgreen3') 
C<-hist(roc1_$Day30_AUC, prob = T) #直方图
lines(density(na.omit(roc1_$Day30_AUC)), col = 'steelblue2') 

load(file = "LASSO结果/seed123-alpha0.33/train55-valid45/ROC2.Rdata")
library(tidyr)  
library(dplyr)  
# 使用list存储重复的AUC值  
df_wide_list <- pivot_wider(roc2, names_from = "time", values_from = "AUC", values_fn = list)  
# 使用平均值汇总重复的AUC值  
df_wide_mean <- pivot_wider(roc2, names_from = "time", values_from = "AUC", values_fn = mean)  
# 识别重复的行  
duplicates <- roc2 %>%  
  group_by(proportion, time) %>%  
  summarise(n = n(), .groups = 'drop') %>%  
  filter(n > 1L)
roc2_ <-df_wide_mean 
library(ggplot2)
library(ggpubr)
library(ggpmisc)
colnames(roc2_)<-c( "Percent",'Order',"Day10_AUC","Day20_AUC","Day30_AUC")
roc2_<-as.data.frame(lapply(roc2_,as.numeric))
colors()
roc2$time<-as.character(roc2$time)
ggplot(aes(proportion,AUC),data = roc2)+
  geom_point(aes(color = time),size = 2.5)+
  scale_x_continuous(limits = c(0.45, 1), breaks = seq(0.45, 1, by = 0.05), labels = scales::percent)+
  scale_color_manual(values=c( "indianred2","mediumaquamarine","skyblue1")) +
  geom_smooth(aes(color = time),method = "gam",formula = y ~ x + I(x ^ 2))+theme_test()+
  theme(axis.title.x = element_text(size = 16,face = "bold",family = "serif"),
        axis.title.y = element_text(size = 16,face = "bold",family = "serif"),
        axis.text.x =  element_text(size = 12,family = "serif"),
        axis.text.y = element_text(size = 12,family = "serif"),
        legend.title = element_text(size = 16,family = "serif"),
        legend.text = element_text(size = 14,family = "serif")) +
  labs(x = "10/20/30 Days AUC (Validation set from 45% to 100%)",
       y = "AUC value",
       color = "Focal species identity")+
  stat_cor(aes(color =time),size = 5.0,label.x = 0.7,label.y = 0.3)+
  theme(legend.position = "none")+
  facet_wrap(~time)
fit <- lm(formula = AUC~proportion, data = roc2)
#计算置信区间，默认 95% 置信区间
confidence <- cbind(roc2, predict(fit, roc2, interval = 'confidence', level = 0.95))
head(confidence)

#计算预测区间，默认 95% 预测区间
prediction <- cbind(roc2, predict(fit, roc2, interval = 'prediction', level = 0.95))
head(prediction)
ggplot(roc2, aes(x = proportion ,y =AUC)) +
  geom_point() +
  scale_x_continuous(limits = c(0.45, 1), breaks = seq(0.45, 1, by = 0.05), labels = scales::percent)+
  geom_ribbon(data = prediction, aes(x = proportion, ymin = lwr, ymax = upr), fill = 'gray', alpha = 0.3) +  #添加上述获得的 95% 预测区间的阴影
  geom_line(data = prediction, aes(y = lwr), color = 'red', linetype = 2) +  #添加上述获得的 95% 预测区间下边界虚线
  geom_line(data = prediction, aes(y = upr), color = 'red', linetype = 2) +  #添加上述获得的 95% 预测区间上边界虚线
  geom_ribbon(data = confidence, aes(ymin = lwr, ymax = upr), fill = 'blue', alpha = 0.3) +  #添加上述获得的  95% 置信区间的阴影
  stat_poly_eq(aes(color =time,label = paste(..eq.label.., ..rr.label.., stat(p.value.label), sep = '~`,`~')), formula = y~poly(x, 1), 
               parse = TRUE,  size = 4.0,label.x = 0.2,label.y = 0.3) +  #添加回归式、R2
  labs(x = "10/20/30 Days AUC (Validation set from 45% to 100%)",
       y = "AUC value")+
  geom_line(data = confidence, aes(y = fit), color = "blue") +  #添加回归线（上述获得的预测均值）
  theme(legend.position = "none")+
  facet_wrap(~time)
