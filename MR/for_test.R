rm (list=ls ())
.libPaths("D:/PackagesR/R_LIBS")

library(TwoSampleMR) 
ao <- available_outcomes()
#新冠相关GWAS
data1<-ao[grepl("COVID-19", ao$trait), ]
data1<-data1[data1$population =="European",]
id1<-c(data1$id)

data<-ao[grepl("Ubiquitin", ao$trait), ]
id<-c(data$id)

#函数
mmrr<-function(x,y){
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
  outcome_dat<-extract_outcome_data(snps=exposure_dat$SNP,
                                    outcomes=y, 
                                    proxies = FALSE,
                                    maf_threshold = 0.01,
                                    access_token = NULL)
  dat <- harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat = outcome_dat
  )
  res <- mr(dat)
  return(res)
  
}
new_mmrr<-function(x,y){
  exposure_datx <-exposure_dat[exposure_dat$id.exposure==x,]
  outcome_daty<-outcome_dat[outcome_dat$id.outcome==y,]
  dat <- harmonise_data(
    exposure_dat = exposure_datx,
    outcome_dat = outcome_daty
  )
  res <- mr(dat)
  return(res)
}

#差距非常小
ress1<-mmrr("ebi-a-GCST010775","ebi-a-GCST90010357" )
#nsnp、P
ress2<-mmrr("ebi-a-GCST010783","prot-a-3144")  
ress3<-mmrr("ebi-a-GCST011076","ebi-a-GCST90010357")
ress4<-mmrr("ebi-a-GCST011078","prot-a-1573")
kk<-new_mmrr("ebi-a-GCST010783","prot-a-3144")

x="ebi-a-GCST010783"
y="prot-a-3144"


newtest<-data.frame()
for (i in id1) {
  for (j in id) {
    print(c(i,j))
    #res<-new_mmrr(i,j)
    res<-mmrr(i,j)
    newtest<-rbind(newtest,res)
  }
}
i<-id1[26]
for (j in id[25:29]) {
     print(c(i,j))
     #res<-new_mmrr(i,j)
       res<-mmrr(i,j)
       newtest<-rbind(newtest,res)
     }

j<-"prot-c-2846_24_2"
save(newtest,file="newtest.Rdata")

newtest_ivw<-newtest[newtest$method=="Inverse variance weighted",]
newtest_iwwp<-newtest_ivw[newtest_ivw$pval<0.05,]
save(newtest_iwwp,file = "need_newtest.Rdata")

load(file = "cov-ubi.Rdata")
load(file = "need2.Rdata")
x<-merge(need2,newtest_iwwp,by="id.exposure")
y<-merge(need2,newtest_iwwp,by="id.outcome")
#id=7、16、23没跑（snp为0
#"ebi-a-GCST010778" "prot-c-2846_24_2"，8,4
#"ebi-a-GCST010778" "prot-c-4474_19_2"，8,27
#"ebi-a-GCST010775" "prot-c-2846_24_2"，18,4，没有snp
#"ebi-a-GCST010775" "prot-c-4474_19_2"，18，27，同
#"ebi-a-GCST010777" "prot-c-2846_24_2"，20,4，同
#"ebi-a-GCST010777" "prot-c-4474_19_2"，20,27，同