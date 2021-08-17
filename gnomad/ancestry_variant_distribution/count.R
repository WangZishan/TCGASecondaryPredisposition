rm(list = ls()); setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition')
tmp<-lapply(c('reshape2','ggplot2','ggrepel'),library,character.only = TRUE)
outdir<-'analysis/data/gnomad/ancestry_variant_distribution/count'

### Load in related information
# The number of patients of different ancestries
tmpcs_count<-read.table('analysis/data/ancestry_info/sta_ancestry_patient.txt',header = T,sep = '\t',stringsAsFactors = F)
cs_count<-tmpcs_count$patient_count; names(cs_count)<-tmpcs_count$ancestry
# variant info from gnomad database
gnomad_variant_ac<-read.table('analysis/data/gnomad/bcftools_process/process/gnomad_variant_info.txt',sep = '\t',stringsAsFactors = F,header = T)
rownames(gnomad_variant_ac)<-gnomad_variant_ac$HGVSg1
# sample-variant info from TCGA
tcga_variant_info<-read.table('analysis/data/variant_data/variants_pca_frq_nonCancer_pathogenic.txt',header = T,sep = '\t',stringsAsFactors = F,quote = '')
tcga_variant_ac<-table(tcga_variant_info$HGVSg1,tcga_variant_info$ancestry)

### Plot
variant_inter<-intersect(rownames(tcga_variant_ac),rownames(gnomad_variant_ac))
(ancestries_inter<-intersect(colnames(tcga_variant_ac),unlist(strsplit(colnames(gnomad_variant_ac),'_'))))
(ancestries_inter<-c('Afr','Amr','Eur','Eas','Sas','Oth'))
Variant_type_color<-c('frameshift'='#66C2A5','missense'='#FC8D62','splice-variant'='#8DA0CB','splice-acceptor'='#A6D854','splice-donor'='#33A02C',
                      'start-lost'='#E78AC3','stop-gain'='#A6761D','stop-lost'='#E5C494','synonymous'='#FFD92F','inframe'='#3A56CE')
for (i in 1:length(ancestries_inter)) {
  tmp<-data.frame(tcga_variant_ac[variant_inter,ancestries_inter[i]],gnomad_variant_ac[variant_inter,paste(ancestries_inter[i],'non_cancer_AC',sep = '_')],stringsAsFactors = F)
  tmp1<-tmp[tmp[,1]!=0 & tmp[,2]!=0,]
  if(nrow(tmp1)!=0){
    tmptcga_variant_ac<-tcga_variant_ac[,setdiff(ancestries_inter,ancestries_inter[i])]
    tmpvariant_other<-rownames(tmptcga_variant_ac)[apply(tmptcga_variant_ac,1,sum)!=0]
    tmpancestry_specific<-rep('Yes',nrow(tmp1)); tmpancestry_specific[rownames(tmp1)%in%tmpvariant_other]<-'No'
    
    tmp2<-data.frame(tmp1,rownames(tmp1),tmpancestry_specific,tcga_variant_info[match(rownames(tmp1),tcga_variant_info$HGVSg1),c('HGVSp','Function','ACMG','Symbol_Merge')],stringsAsFactors = F)
    colnames(tmp2)[1:4]<-c('TCGA_AC','Gnomad_AC','HGVSg1','Specificity')
    tmp2$Label<-paste(tmp2$Symbol_Merge,sapply(strsplit(tmp2$HGVSp,':'),function(x){x[2]}),sep = ':')
    
    tmpfun<-function(x,n){if(length(x)>n){return(which(x>=sort(x,decreasing = T)[n]))}else{return(1:length(x))}}#for function
    if(sum(tmp2$Specificity=='Yes')!=0){
      tmptmp2<-tmp2[tmp2$Specificity=='Yes',]
      tmptmp3<-unique(rbind(tmptmp2[union(tmpfun(tmptmp2$TCGA_AC,1),tmpfun(tmptmp2$Gnomad_AC,1)),],tmp2[union(tmpfun(tmp2$TCGA_AC,1),tmpfun(tmp2$Gnomad_AC,1)),]))
    }else{
      tmptmp3<-unique(tmp2[union(tmpfun(tmp2$TCGA_AC,1),tmpfun(tmp2$Gnomad_AC,1)),])
    }#for if

    p<-ggplot(tmp2,aes(x=TCGA_AC,y=Gnomad_AC,color=Function,shape=Specificity))+geom_point(size=3,alpha=0.4)
    p<-p+geom_label_repel(aes(label=ifelse(HGVSg1 %in% rownames(tmptmp3), Label, NA)),size=5)
    p<-p+scale_color_manual(values = Variant_type_color)
    p<-p+labs(title=paste(ancestries_inter[i],'AC'))
    p<-p+theme_bw()
    p<-p+theme(axis.title = element_text(size = 22),axis.text.y = element_text(size = 20),
               axis.text.x = element_text(size = 20),
               legend.title = element_text(size = 22),legend.text = element_text(size = 20),
               legend.key.size = unit(1.3,'cm'),legend.box.background = element_blank(),
               panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),
               title = element_text(size = 25))
    ylim_val<-max(tmp2$Gnomad_AC)
    if(max(tmp2$TCGA_AC)<(ylim_val/3)){xlim_val=ylim_val/3}else{xlim_val<-max(tmp2$TCGA_AC)}#for if
    p<-p+coord_fixed(ratio = 1,xlim = c(0,xlim_val*1.2),ylim = c(0,ylim_val*1.2))
    ggsave(p,file=file.path(outdir,paste(ancestries_inter[i],'_AC.pdf',sep = '')),h=8,w=10,useDingbat=F)
    
    tmp2$Labelled<-rownames(tmp2)%in%rownames(tmptmp3)
    if(nrow(tmp2)==1){tmptmp31<-t(tcga_variant_ac[rownames(tmp2),ancestries_inter])}else{tmptmp31<-tcga_variant_ac[rownames(tmp2),ancestries_inter]}#for if
    tmptmp32<-t(t(tmptmp31)/cs_count[ancestries_inter])
    tmptmp33<-cbind(tmptmp32,tmptmp31)
    colnames(tmptmp33)<-c(paste('TCGA_',ancestries_inter,'_AF',sep = ''),paste('TCGA_',ancestries_inter,'_AC',sep = ''))
    tmp3<-cbind(tmp2[,c('Specificity','Labelled','Symbol_Merge','HGVSp','Function','ACMG')],tmptmp33,
                gnomad_variant_ac[rownames(tmp2),c(paste(ancestries_inter,'non_cancer_AF',sep = '_'),paste(ancestries_inter,'non_cancer_AC',sep = '_'))])
    tmp4<-tmp3[order(tmp3$Specificity,tmp3[,paste('TCGA_',ancestries_inter[i],'_AC',sep = '')],decreasing = T),]
    write.table(tmp4,file.path(outdir,paste(ancestries_inter[i],'_AC.txt',sep = '')),row.names = T,col.names = T,quote = F,sep = '\t')
  }#for if
}#for i

