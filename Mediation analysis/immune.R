
rm(list=ls())

##### 中介孟德尔随机化 #####

library(TwoSampleMR)
library(kableExtra)
library(ggplot2)
library(cowplot)
library(dplyr) 
library(tidyr)
library(RMediation)

load(file = "ao.Rdata")
ao<-as.data.frame(ao)
X<-"ebi-a-GCST011078"
Y<-"prot-a-2592"

#X->Y 正向
exposure_dat1 <-extract_instruments(outcomes=X,p1=5e-6,clump=TRUE, r2=0.001,kb=10000)
Ffilter=10        #F值过滤条件
d<-exposure_dat1 
#计算F检验值
N=d[1,"samplesize.exposure"]     #获取样品的数目
d=transform(d,R2=2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))     #计算R2
d=transform(d,F=(N-2)*R2/(1-R2))      #计算F检验值
#根据F值>10进行过滤, 删除弱工具变量
outTab=d[d$F>Ffilter,]
cols <- ncol(outTab)
is_na <- is.na(outTab)
row_na <- rowSums(is_na)
outTab<- outTab[row_na != cols, ]
exposure_dat1<-outTab
outcome_dat<-extract_outcome_data(snps=exposure_dat1$SNP,
                                  outcomes=Y, 
                                  proxies = FALSE,
                                  maf_threshold = 0.01)
dat <- harmonise_data(
  exposure_dat = exposure_dat1,
  outcome_dat = outcome_dat
)
res <- mr(dat)

#Y->X （反向）
exposure_dat2 <-extract_instruments(outcomes=Y,p1=5e-6,clump=TRUE, r2=0.001,kb=10000)
Ffilter=10        #F值过滤条件
d<-exposure_dat2 
#计算F检验值
N=d[1,"samplesize.exposure"]     #获取样品的数目
d=transform(d,R2=2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))     #计算R2
d=transform(d,F=(N-2)*R2/(1-R2))      #计算F检验值
#根据F值>10进行过滤, 删除弱工具变量
outTab=d[d$F>Ffilter,]
cols <- ncol(outTab)
is_na <- is.na(outTab)
row_na <- rowSums(is_na)
outTab<- outTab[row_na != cols, ]
exposure_dat2<-outTab
outcome_dat<-extract_outcome_data(snps=exposure_dat2$SNP,
                                  outcomes=X, 
                                  proxies = FALSE,
                                  maf_threshold = 0.01)
dat <- harmonise_data(
  exposure_dat = exposure_dat2,
  outcome_dat = outcome_dat
)
res_reverse <- mr(dat)

#  1) Delta
product_method_Delta <- function(EM_beta, EM_se, MO_beta, MO_se, verbose=F){
 
  # 计算中介效应beta
  EO <- EM_beta * MO_beta
  
  if (verbose) {print(paste("Indirect effect = ", round(EM_beta, 2)," x ", round(MO_beta,2), " = ", round(EO, 3)))}
  
  
  #  RMediation package计算置信区间
  CIs = medci(EM_beta, MO_beta, EM_se, MO_se, type="dop")
  
  # 计算 Z 值和 p 值
  z_score <- EO / CIs$SE
  p_value <- 2 * pnorm(abs(z_score), lower.tail = FALSE)
  
  # 把数据 放进数据集
  df <-data.frame(Mediated_effect = EO,
                  #Mediated_se = CIs$SE,
                  Mediated_p = p_value,
                  CI = paste("[",round(CIs[["95% CI"]][1],4),",",
                             round(CIs[["95% CI"]][2],4),"]"))

  return(df)
}

#"beta_total","p_total","beta_reverse","p_reverse",                            
#"Mediated effect","CI","p","Mediated proption"
immid <- paste0('ebi-a-GCST9000',c(1391:2121))
mediation_T<-data.frame()
mediation_F<-data.frame()
for (M in immid[(sum(nrow(mediation_T)+nrow(mediation_F))+1):731]) {
  test<-data.frame("M"=M)
  test$trait<-ao$trait[ao$id == M]
  #X->M
  outcome_dat1 <- NULL 
  while(is.null(outcome_dat1) || nrow(outcome_dat1) == 0) {
    try({
      outcome_dat1<-extract_outcome_data(snps=exposure_dat1$SNP,
                                         outcomes=M, 
                                         proxies = FALSE,
                                         maf_threshold = 0.01)
    })
    Sys.sleep(2)
  }
  
  dat1 <- harmonise_data(
    exposure_dat = exposure_dat1,
    outcome_dat = outcome_dat1
  )
  res1 <- mr(dat1)
  
  test$beta_EM<-res1 %>%filter(method == "Inverse variance weighted") %>%pull(b)
  test$se_EM<-res1 %>%filter(method == "Inverse variance weighted")%>% pull(se)
  test$p_EM<-res1 %>%filter(method == "Inverse variance weighted")%>% pull(pval)
  
  #M
  exposure_dat3 <- NULL 
  while(is.null(exposure_dat3) || nrow(exposure_dat3) == 0) {
    try({
      exposure_dat3 <-extract_instruments(outcomes=M,p1=5e-6,clump=TRUE, r2=0.001,kb=10000)
    })
    Sys.sleep(2)
  }
  R<- TwoSampleMR::get_r_from_bsen(exposure_dat3$beta.exposure,exposure_dat3$se.exposure,exposure_dat3$samplesize.exposure)
  if (!is.na(R)[1]) {
    F_value <- as.numeric(R*R)*(as.numeric(exposure_dat3$samplesize.exposure) - 2)/(1 - as.numeric(R*R))
    exposure_dat3<- cbind(exposure_dat3, R2 = as.numeric(R*R), F_value = F_value)
    exposure_dat3 <- exposure_dat3[as.numeric(exposure_dat3$F_value) > 10,]
  } 
  #M->Y
  outcome_dat2<-NULL 
  while(is.null(outcome_dat2) || nrow(outcome_dat2) == 0) {
    try({
      outcome_dat2<-extract_outcome_data(snps=exposure_dat3$SNP,
                                         outcomes=Y, 
                                         proxies = FALSE,
                                         maf_threshold = 0.01)
    })
    Sys.sleep(2)
  }
  dat2 <- harmonise_data(
    exposure_dat = exposure_dat3,
    outcome_dat = outcome_dat2
  )
  res2 <- mr(dat2)
  tryCatch({  
    test$beta_MO<-res2 %>%filter(method == "Inverse variance weighted") %>%pull(b)
    test$se_MO<-res2 %>%filter(method == "Inverse variance weighted")%>% pull(se)
    test$p_MO<-res2 %>%filter(method == "Inverse variance weighted")%>% pull(pval)
  }, error = function(e) { 
    test$beta_MO <<- NA
    test$se_MO <<- NA 
    test$p_MO <<- NA })  
  if(test$p_EM<0.05 && test$p_MO<0.05){
    Delta<-product_method_Delta(test$beta_EM, test$se_EM, test$beta_MO, test$se_MO)
    test<-cbind(test,Delta) 
    test$Mediated_proption<-test$Mediated_effect/res[res$method=="Inverse variance weighted",]$b
    mediation_T<-rbind(mediation_T,test)
  }else
    mediation_F<-rbind(mediation_F,test)
}

save(mediation_F,mediation_T,file = "mediation_immune.Rdata")

Total_effect<-res %>%filter(method == "Inverse variance weighted") %>%pull(b)
mediation_T$Direct_effect<-Total_effect-mediation_T$Mediated_effect

#save(mediation_F,mediation_T,file = "mediation_immune_e-5.Rdata")
load(file = "mediation_immune_5_10e-6.Rdata")
write.csv(mediation_F,file = "mediation_F.csv")
write.csv(mediation_T,file = "mediation_T.csv")
