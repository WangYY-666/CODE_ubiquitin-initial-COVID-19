#整体计算
rm (list=ls ())
.libPaths("D:/PackagesR/R_LIBS")

library(TwoSampleMR) 
load(file = "need_newtest.Rdata")
x=newtest_iwwp$id.exposure[1]
y=newtest_iwwp$id.outcome[1]
#暴露
exposure_dat <-extract_instruments(outcomes=x,p1=5e-6,clump=TRUE, r2=0.001,kb=10000,access_token= NULL)

Ffilter=10        #F值过滤条件
d<-exposure_dat 
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
exposure_dat<-outTab

#结局
outcome_dat<-extract_outcome_data(snps=exposure_dat$SNP,
                                  outcomes="prot-a-3181", 
                                  proxies = FALSE,
                                  maf_threshold = 0.01,
                                  access_token = NULL)
#整合
dat <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat = outcome_dat
)
#save(dat,file = "dat.Rdata")
load(file = "dat.Rdata")
#孟德尔随机化
res <- mr(dat)

#计算系数的置信区间及OR值
res_or<-generate_odds_ratios(res)
save(res_or,file = "res(and OR).Rdata")

library(writexl)
#异质性检验
a<-mr_heterogeneity(dat)
write_xlsx(a,"异质性检验.xlsx")


#水平多效性检验
#如果变量工具不通过暴露影响结果，就违反了孟德尔的假设，就是存在多水平效应。
b<-mr_pleiotropy_test(dat)
write_xlsx(b,"水平多效性.xlsx")

#MR Steiger 方向性测试
#在 MR 中，假定工具首先影响暴露，然后通过暴露影响结果。但有时这很难评估，例如顺式作用 SNP 是先影响基因表达水平还是先影响 DNA 甲基化水平？这时就可以使用 Steiger 检验来检验假设暴露与结果之间的因果方向。
out <- directionality_test(dat)

x11()
#留一法（一种交叉验证方法）
res_loo <- mr_leaveoneout(dat)
pdf("留一法敏感性检验.pdf",width=6,height=7.5)
mr_leaveoneout_plot(res_loo)
dev.off()

#散点图
pdf("散点图.pdf")
p1 <- mr_scatter_plot(res, dat)
p1
dev.off()

#森林图
pdf("森林图.pdf")
res_single <- mr_singlesnp(dat)
mr_forest_plot(res_single)
dev.off()

#漏斗图
pdf("漏斗图.pdf")
c<-mr_funnel_plot(res_single)
c
dev.off()

#直接生成报告单
mr_report(dat)
