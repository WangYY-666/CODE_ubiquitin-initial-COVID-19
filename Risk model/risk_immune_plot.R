
##续IOBR

res<-cell_bar_plot(input = im_cibersort[1:12,], features = colnames(im_cibersort)[3:22], title = "CIBERSORT Cell Fraction")

##
cols <- IOBR::palettes(category = "random", palette = 4,
                       show_col = F,
                       show_message = F)
plot_tme = function(meltdf){
  p_stack_df = ggplot(meltdf,aes(x = sample, y = value, fill = variable)) + 
    geom_bar(stat = "identity")+ theme_test() + 
    scale_fill_manual(values = cols) + #scale_x_discrete(limits = rev(levels(input))) + 
    #ggtitle(paste0(title)) + 
    theme(legend.position = 'bottom')+
    theme(plot.title = element_text(size = rel(2),hjust = 0.5), 
          axis.text.x = element_blank())+
    theme(axis.ticks =element_blank() )+
    scale_y_continuous(expand = c(0,0))+
    facet_nested(.~status1,drop=T,scale="free",space="free",switch="x",
                 strip =strip_nested(background_x = elem_list_rect(fill =c('steelblue','firebrick')),by_layer_x = F))
  
  # p_stack_df
  
  p_box_df = ggplot(meltdf,aes(x = variable, y = value, 
                               fill = variable)) +
    geom_boxplot(width=0.4,lwd=0.2,color='black',outlier.shape=NA,
                 position = position_dodge(width = 0.8))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))+
    theme(legend.position = 'none')+
    scale_fill_manual(values = cols)
  # p_box_df
  
  p_box_comp =  ggplot(meltdf,aes(x = variable, y = value, 
                                  fill = status1,color = status1)) +
    
    geom_boxplot(width=0.5,lwd=0.2,color='black',outlier.shape=NA,
                 position = position_dodge(width = 0.8))+
    theme_bw()+
    theme(legend.position = 'top')+
    theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))+
    stat_compare_means(aes(group=status1,label = after_stat(p.signif)),hide.ns = FALSE,
                       method = "t.test")+
    scale_fill_manual(values = c('steelblue','firebrick'))
  
  # p_box_comp
  plot_comb = (p_stack_df)/(p_box_df|p_box_comp)
  return(plot_comb)

}
####cibersort
raw_ciber<-im_cibersort
raw_ciber<-  column_to_rownames(raw_ciber, "ID")
names(raw_ciber)<-sub("_CIBERSORT$", "", names(raw_ciber))  
raw_ciber<-cbind(ID=rownames(raw_ciber),raw_ciber)
re_ciber= raw_ciber%>%
  as.data.frame()%>%
  left_join(cli%>%dplyr::select(c('ID','risk')),by='ID')%>%
  dplyr::select(c('ID','risk'),everything())%>%
  #dplyr::select(-c((ncol(.) - 2):ncol(.)))%>%
  arrange((risk))

meltdf_ciber = re_ciber%>%
  reshape2::melt()%>%
  mutate(ID=factor(ID,levels=re_ciber$ID))%>%
  mutate(risk=factor(risk,levels=c('low','high')))
colnames(meltdf_ciber)<-c("sample","status1","variable", "value" )
library(ggh4x)
plot_ciber= plot_tme(meltdf_ciber)
print(plot_ciber)
pdf('./plot_comb_ciber1.pdf',width = 14,height = 14)
dev.off()
####ssgsea
raw<-im_ssgsea
re_ssgsea= raw%>%
  as.data.frame()%>%
  left_join(cli%>%dplyr::select(c('ID','risk')),by='ID')%>%
  dplyr::select(c('ID','risk'),everything())%>%
  #dplyr::select(-c((ncol(.) - 2):ncol(.)))%>%
  arrange((risk))

meltdf= re_ssgsea%>%
  reshape2::melt()%>%
  mutate(ID=factor(ID,levels=re_ciber$ID))%>%
  mutate(risk=factor(risk,levels=c('low','high')))
colnames(meltdf)<-c("sample","status1","variable", "value" )
library(ggh4x)
plot_ssgsea= plot_tme(meltdf)
print(plot_ssgsea)


###xcell
res_xcell<-  column_to_rownames(im_xcell, "ID")
names(res_xcell)<-sub("_xCell$", "", names(res_xcell))  
res_xcell<-cbind(ID=rownames(res_xcell),res_xcell)

load(file = "risksocre_order.Rdata")
cli<-rt %>%
  select("order" ,"risk")
colnames(cli)[1]<-"ID"

res_bind<-merge(cli,res_xcell,by="ID")
rownames(res_bind)<-res_bind[,1]
res_bind<-res_bind[,-1]

colnames(res_bind)[1]="Group"

rt<-res_bind
library(tidyverse)

rt2=t(rt[,-1])
group=rt[,1]
phylum=as.data.frame(rt2)
phylum_per <- as.data.frame(lapply(phylum, function(x) x / sum(x)))
row.names(phylum_per) <- row.names(phylum)

phylum.ave <- apply(phylum_per, 1, FUN=mean)
phylum.2 <- cbind(phylum_per, phylum.ave)[order(-phylum.ave),]


phylum.2 <- subset(phylum.2[1:10,], select=-phylum.ave)
phylum.2 <- rbind(phylum.2, others=apply(phylum.2, 2, function(x){1-sum(x)}))
library(reshape2)

PhylumID=row.names(phylum.2)
phylum.2=cbind(PhylumID,phylum.2)
phylum.2$PhylumID <- factor(phylum.2$PhylumID, levels = rev(phylum.2$PhylumID))
phylum.gg <- melt(phylum.2, id.vars="PhylumID", variable.name="SampleID", value.name="Abundance")
head(phylum.gg)
phylum.gg$group <- group



library(wesanderson)
library(colortools)
library(ggpubr)

ggbarplot(phylum.gg, x = "SampleID", y="Abundance", color="black", fill="PhylumID",
          legend="right", 
          legend.title="Top immune cells", main="Relative immune cells per patient by xCell",
          font.main = c(12,"bold", "black"), font.x = c(12, "bold"), 
          font.y=c(12,"bold")) + 
  theme_classic()  +
  rotate_x_text() + 
  scale_fill_manual(values=c("gray",rev(wheel("skyblue3")[1:10]))) +
  facet_grid(~ group, scales = "free_x", space='free',) + 
  labs(x = "Sample", y = "Relative proportion") + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_text(size = 20,vjust = 0.5,hjust = 0.5,face = "bold"), 
        plot.title = element_text(size = 20,vjust = 0.5,hjust = 0.5,face = "bold"), 
        legend.title = element_text(size = 20,vjust = 0.5,hjust = 0.5,face = "bold"),
        axis.title.y =  element_text(size = 20,vjust = 0.5,hjust = 0.5,face = "bold"),
        legend.text =  element_text(size = 10,face = "bold")
  ) 
ggsave(filename = "xCell_topimmune.pdf", device="pdf", width=20, height=6)

