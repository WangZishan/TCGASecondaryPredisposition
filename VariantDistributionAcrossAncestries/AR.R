rm(list = ls()); tmp<-lapply(c('ggplot2','reshape2','grid'), library,character.only=T)
setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition/')

### load in related information
# variant info
variant<-read.table('analysis/data/variant_data/variants_pca_frq_nonCancer_pathogenic.txt',header = T,sep = '\t',quote = '',stringsAsFactors = F)
variant <- variant[variant$Chromosome!='X',]
variant0<-unique(variant[,c('bcr_patient_barcode','ancestry','Symbol_Merge','AR')])
tmpind<-grep(';',variant0$Symbol_Merge)
if(length(tmpind)!=0){
  tmpvariant1<-variant0[setdiff(1:nrow(variant0),tmpind),]; tmpvariant2<-variant0[tmpind,]
  tmpvariant2_dup<-matrix(nrow = 0,ncol = ncol(variant0))
  for (i in 1:nrow(tmpvariant2)) {
    tmp0<-matrix(nrow = 0,ncol = ncol(variant0)); tmp1<-unlist(strsplit(tmpvariant2[i,'Symbol_Merge'],';'))
    for (j in 1:length(tmp1)) {tmp0<-rbind(tmp0,tmpvariant2[i,])}# for
    tmp0$Symbol_Merge<-tmp1; tmpvariant2_dup<-rbind(tmpvariant2_dup,tmp0)
  }# for
  variant1<-unique(rbind(tmpvariant1,tmpvariant2_dup))
}else{
  variant1 <- unique(variant0)
}
# ancestry info
ancestry<-read.table('analysis/data/ancestry_info/ancestry_patient.txt',header = T,sep = '\t',stringsAsFactors = F)
ancestry_sta<-table(ancestry$ancestry); ancestry_order<-c('Amr','Afr','Eur','Eas','Sas','Mix','Oth')

###count/proportion of patients among AR status across ancestries
# count/proportion statistic
tmpvariant<-expand.grid(ancestry=unique(variant$ancestry),AR=unique(variant$AR),Count1=NA,Count2=NA,Proportion1=NA,Proportion2=NA,stringsAsFactors = F)
for (i in 1:nrow(tmpvariant)) {
  tmpvariant[i,'Count1']<-sum(variant$ancestry==tmpvariant[i,'ancestry'] & variant$AR==tmpvariant[i,'AR'])
  tmpvariant[i,'Count2']<-sum(variant$ancestry==tmpvariant[i,'ancestry'])
  if(tmpvariant[i,'AR']=='No'){
    tmpvariant[i,'Proportion1']<-length(unique(variant[(variant$ancestry==tmpvariant[i,'ancestry'] & variant$AR=='No'),'bcr_patient_barcode']))/ancestry_sta[tmpvariant[i,'ancestry']]
  }else if(tmpvariant[i,'AR']=='Yes'){
    tmp <- unique(setdiff(variant[(variant$ancestry==tmpvariant[i,'ancestry'] & variant$AR=='Yes'),'bcr_patient_barcode'],variant[(variant$ancestry==tmpvariant[i,'ancestry'] & variant$AR=='No'),'bcr_patient_barcode']))
    tmpvariant[i,'Proportion1']<-length(tmp)/ancestry_sta[tmpvariant[i,'ancestry']]
  }
  tmpvariant[i,'Proportion2']<-length(unique(variant[(variant$ancestry==tmpvariant[i,'ancestry']),'bcr_patient_barcode']))/ancestry_sta[tmpvariant[i,'ancestry']]
}# for
write.table(tmpvariant,'analysis/variant_distribution_of_ancestry/AR_distribution/AR.txt',row.names = F,col.names = T,quote = F,sep = '\t')
tmpvariant$ancestry<-factor(tmpvariant$ancestry,levels = rev(ancestry_order))
tmpvariant1<-tmpvariant[tmpvariant$AR=='Yes',];tmpvariant2<-tmpvariant[tmpvariant$AR=='No',]
# count plot
p<-ggplot()+geom_bar(data = tmpvariant2, aes(x=ancestry,y=Count2),stat = 'identity',fill='#F8766D')
p<-p+geom_bar(data = tmpvariant1, aes(x=ancestry,y=Count1),stat = 'identity',fill='#00BFC4')
p<-p+geom_text(data = tmpvariant2,aes(x=ancestry,y=Count2*0.65,label=Count1),size=6)
p<-p+geom_text(data = tmpvariant1,aes(x=ancestry,y=1.8,label=Count1),size=6)
p<-p+theme_bw()+ theme(axis.text = element_text(size = 20),axis.ticks.length.y = unit(0,'cm'),
                                          axis.title.x =  element_text(size = 23),axis.title.y = element_blank(),
                                          legend.title = element_text(size = 23),legend.text = element_text(size = 20),
                                          panel.grid.major.y = element_blank(),legend.position = 'none')
p<-p+labs(x='Ancestry',y='Secondary findings count')+coord_flip()+scale_y_continuous(trans = 'log10')
p
ggsave('analysis/variant_distribution_of_ancestry/AR_distribution/AR_Count.pdf',h=6,w=3,useDingbat=F)
#proportion plot
p<-ggplot()+geom_bar(data = tmpvariant2, aes(x=ancestry,y=Proportion2),stat = 'identity',fill='#F8766D')
p<-p+geom_bar(data = tmpvariant1, aes(x=ancestry,y=Proportion1),stat = 'identity',fill='#00BFC4')
p<-p+geom_text(data = tmpvariant2,aes(x=ancestry,y=Proportion2+0.05,label=round(tmpvariant2$Proportion1*100,2)),size=6)
p<-p+geom_text(data = tmpvariant1,aes(x=ancestry,y=0.05,label=round(tmpvariant1$Proportion1*100,2)),size=6)
p<-p+labs(x='Ancestry',y='Carrier Proportion(%)')+coord_flip()+scale_y_reverse()
p<-p+theme_bw()+ theme(axis.text = element_text(size = 20),axis.ticks.length.y = unit(0,'cm'),
                       axis.title.x =  element_text(size = 23),axis.title.y = element_blank(),
                       legend.title = element_text(size = 23),legend.text = element_text(size = 20),
                       panel.grid.major.y = element_blank(),legend.position = 'none')
p
ggsave('analysis/variant_distribution_of_ancestry/AR_distribution/AR_Proportion.pdf',h=6,w=3,useDingbat=F)


