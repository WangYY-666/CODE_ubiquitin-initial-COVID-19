rm(list = ls())
.libPaths("D:/PackagesR/R_LIBS")

load(file = "GSE157103_survivaltest.Rdata")
library("survival")
library("survminer")
colnames(rt_cli)[2:3]<-c("mechanical_ventilation","ventilator_free_days")
expr<-as.data.frame(lapply(rt_cli,as.numeric))
expr_normal<-expr[,4:19475]
expr_normal<-log2(expr_normal+1)
expr[,4:19475]<-expr_normal
expr[,1]<-rt_cli[,1]
save(expr,file = "GSE157103_survivaltest_normalized.Rdata")
#单因素cox回归分析
#res.cox <- coxph(Surv(ventilator_free_days,mechanical_ventilation) ~ A1CF, data =expr)
#res.cox
#这里的exp(coef)就是HR（hazard ratio，风险率）
#lower .95和upper .95为95%的置信区间
#summary(res.cox)# summary查看更多信息
##对关注的单基因进行cox回归：
#gene <- c("USP8")
#x<-expr[,c(1,4:19475)]
#转换为二分类变量(高低表达两组)：
#expr$gene<- ifelse(x[,gene]>median(x[,gene]),"High","Low")
#cox_USP8 <- coxph(Surv(ventilator_free_days,mechanical_ventilation) ~ gene, data = expr)
#cox_USP8
#summary(cox_USP8)
###

gene<-data.table::fread("gene_need.txt")
#假设我们要对如下特征做单因素cox回归分析
covariates <- c("USP8" , "USP25","USP21","UFM1" , "UBL4A","UBE2V1","UBE2N", "UBE2L3" ,"UBE2G2" ,"UBE2B" , "OTUB2" , "ATG7", "ISG15" )
#分别对每一个变量，构建生存分析的公式
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(ventilator_free_days,mechanical_ventilation)~', x)))
#循环对每一个特征做cox回归分析
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = expr)})
#提取HR，95%置信区间和p值
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         #获取HR
                         HR <-signif(x$coef[2], digits=2);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })
#转换成数据框，并转置
res <- t(as.data.frame(univ_results, check.names = FALSE))
res<-as.data.frame(res)

res_p<-res[res$p.value<0.05,]
write.table(file="单因素cox结果.txt",as.data.frame(res),quote=F,sep="\t")
save(res_p,file = "res_p.Rdata")
######

#多因素cox回归分析
load(file = "res_p.Rdata")
rownames(res_p)
load(file = "GSE157103_survivaltest_normalized.Rdata")

#一般先通过单因素cox回归分析找出与生存显著相关的特征，然后基于这些特征再去做多因素cox回归分析，或者做LASSO分析。
res.cox <- coxph(Surv(ventilator_free_days,mechanical_ventilation) ~ USP21+UBE2V1+UBE2N+UBE2L3+UBE2G2+UBE2B+ATG7+ISG15, data = expr)
x <- summary(res.cox)
pvalue=signif(as.matrix(x$coefficients)[,5],2)
HR=signif(as.matrix(x$coefficients)[,2],2)
low=signif(x$conf.int[,3],2)
high=signif(x$conf.int[,4],2)
multi_res=data.frame(p.value=pvalue,
                     HR=paste(HR," (",low,"-",high,")",sep=""),
                     stringsAsFactors = F
)
multi_res
write.table(file="多因素cox结果.txt",multi_res,quote=F,sep="\t")
multi_res_p<-multi_res[multi_res$p.value<0.05,]
multi_res_p
#其他算法
library("autoReg")
#cox回归模型构建
coxmod<-coxph(Surv(ventilator_free_days,mechanical_ventilation)~USP8+USP25+USP21+UFM1+UBL4A+UBE2V1+UBE2N+UBE2L3+UBE2G2+UBE2B+OTUB2+ATG7+ISG15,
              data=expr)
summary(coxmod)
#单因素＋P<0.05纳入多因素＋逐步回归后退法
ft3<-autoReg(coxmod,uni=TRUE,threshold=0.05, final= TRUE) 
myft(ft3)



#https://mp.weixin.qq.com/s/0bCT6zBGNBSpalhjwPx-CQ

#森林图
options(scipen=1)
x11()
ggforest(res.cox,
         data = expr,
         main = "Hazard ratio",
         cpositions = c(0.02, 0.22, 0.4),
         fontsize = 1,
         refLabel = "1",
         noDigits = 3)
#C-inddex计算
summary(res.cox)#提取Concordance
C_index <- summary(res.cox)$concordance[1]
C_index
#计算风险评分(是每组的风险评分)
gene<-c("USP21" ,"UBE2V1", "UBE2N" , "UBE2L3" ,"UBE2G2" ,"UBE2B" ,"ATG7","ISG15"  )
riskScore <- predict(res.cox,type="risk",newdata=expr[,gene]) 

#绘制ROC
load(file = "GSE157103_survivaltest_normalized.Rdata")
df<-expr
library(multipleROC)
###多个ROC曲线绘制
p1 <- multipleROC(mechanical_ventilation~UBE2V1,data=df)
p2 <- multipleROC(mechanical_ventilation~UBE2L3,data=df)
p3 <- multipleROC(mechanical_ventilation~UBE2B,data=df)
p4 <- multipleROC(mechanical_ventilation~ISG15,data=df)
#一起显示
plot_ROC(list(p1,p2,p3,p4),
         show.points = T, 
         show.eta = F, 
         show.sens = F, 
         show.AUC = T, 
         facet = F )
#分面显示
plot_ROC(list(p1,p2,p3,p4),
         show.points = T, 
         show.eta = T, 
         show.sens = T, 
         show.AUC = T, 
         facet = T )
library(pROC)
res<-roc(mechanical_ventilation~UBE2V1+UBE2L3+UBE2B+ISG15,data=df,aur=TRUE,
         ci=TRUE, # 显示95%CI
         #percent=TRUE, # 是否需要以百分比显示
         smooth=TRUE,# 是否平滑曲线
)
library(ggplot2)
p<- ggroc(res, legacy.axes = TRUE)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
  theme_bw() + # 设置背景
  ggtitle("ROC Curve")+
  theme(plot.title = element_text(hjust = 0.5,size = 16),
        axis.text=element_text(size=12,colour = "black"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

p+annotate("text",x=0.65,y=0.05,label=paste("UBE2V1-AUC = ", round(res$UBE2V1$auc,3)))+
  annotate("text",x=0.65,y=0.15,label=paste("UBE2L3-AUC = ", round(res$UBE2L3$auc,3)))+
  annotate("text",x=0.65,y=0.25,label=paste("UBE2B-AUC = ", round(res$UBE2B$auc,3)))+
  annotate("text",x=0.65,y=0.35,label=paste("ISG15-AUC = ", round(res$ISG15$auc,3)))

ROC1<-roc(mechanical_ventilation~UBE2V1,data=df,aur=TRUE,
         ci=TRUE, # 显示95%CI
         #percent=TRUE, # 是否需要以百分比显示
         smooth=TRUE,# 是否平滑曲线
)

ROC2<-roc(mechanical_ventilation~UBE2B,data=df,aur=TRUE,
          ci=TRUE, # 显示95%CI
          #percent=TRUE, # 是否需要以百分比显示
          smooth=TRUE,# 是否平滑曲线
)
ROC3<-roc(mechanical_ventilation~UBE2L3,data=df,aur=TRUE,
          ci=TRUE, # 显示95%CI
          #percent=TRUE, # 是否需要以百分比显示
          smooth=TRUE,# 是否平滑曲线
)
ROC4<-roc(mechanical_ventilation~ISG15,data=df,aur=TRUE,
          ci=TRUE, # 显示95%CI
          #percent=TRUE, # 是否需要以百分比显示
          smooth=TRUE,# 是否平滑曲线
)
plot(ROC1, col = "red2", lwd = 2) #绘制ROC曲线，定义曲线颜色和线条粗细
plot(ROC2, col = "blue2", lwd = 2,add = TRUE)  
plot(ROC3, col = "green3", lwd = 2,add = TRUE)  
plot(ROC4, col = "yellow2", lwd = 2,add = TRUE)  
ROC1
ROC2
ROC3
ROC4
legend("bottomright", cex = 1.1,
       legend = c("UBE2V1 (AUC: 0.783)", "UBE2B (AUC: 0.728)", "UBE2L3 (AUC: 0.696)","ISG15 (AUC: 0.776)"),
       col = c("red2", "blue2", "green3","yellow2"), lty = 1, lwd = 2)  #定义标签

#生存分析
rownames(multi_res)
ggsurvplot(survfit(res.cox, data = expr), palette = "#2E9FDF")
