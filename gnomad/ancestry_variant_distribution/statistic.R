rm(list = ls()); setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition')

### load in related information
# ancestry info
ancestry<-read.table('analysis/data/ancestry_info/ancestry_patient.txt',header = T,sep = '\t',stringsAsFactors = F)
# variant info from gnomad database
gnomad_variant<-read.table('analysis/data/gnomad/bcftools_process/process/gnomad_variant_info.txt',sep = '\t',stringsAsFactors = F,header = T)
rownames(gnomad_variant)<- gnomad_variant$HGVSg1
# sample-variant info from TCGA
tcga_variant_info<-read.table('analysis/data/variant_data/variants_pca_frq_nonCancer_pathogenic.txt',header = T,sep = '\t',stringsAsFactors = F,quote = '')
tcga_variant_ac<-table(tcga_variant_info$HGVSg1,tcga_variant_info$ancestry)
tcga_variant_af<-t(t(tcga_variant_ac)/as.numeric(table(ancestry$ancestry)[colnames(tcga_variant_ac)]))

### Calculating information of interests
variant_inter<-intersect(rownames(tcga_variant_ac),rownames(gnomad_variant))
ancestries_inter<-intersect(colnames(tcga_variant_ac),unlist(strsplit(colnames(gnomad_variant),'_')))
m<-matrix(nrow = length(ancestries_inter),ncol = 12,dimnames = list(ancestries_inter,c('Num_plot','Num_spe','Num_not_spe','Num_label',paste('Freq',c('pearson','pearson_p','spearman','spearman_p'),sep = '_'),paste('Count',c('pearson','pearson_p','spearman','spearman_p'),sep = '_'))))
for (i in 1:length(ancestries_inter)) {
  tmp<-data.frame(tcga_variant_ac[variant_inter,ancestries_inter[i]],gnomad_variant[variant_inter,paste(ancestries_inter[i],'non_cancer_AC',sep = '_')],
                  tcga_variant_af[variant_inter,ancestries_inter[i]],gnomad_variant[variant_inter,paste(ancestries_inter[i],'non_cancer_AF',sep = '_')],stringsAsFactors = F)
  tmp1<-tmp[(tmp[,1]!=0) & (tmp[,2]!=0),]
  if(nrow(tmp1)!=0){
    tmptcga_variant_ac<-tcga_variant_ac[,setdiff(ancestries_inter,ancestries_inter[i])]
    tmpvariant_other<-rownames(tmptcga_variant_ac)[apply(tmptcga_variant_ac,1,sum)!=0]
    tmpancestry_specific<-rep('Yes',nrow(tmp1))
    tmpancestry_specific[rownames(tmp1)%in%tmpvariant_other]<-'No'
    
    tmp2<-data.frame(tmp1,rownames(tmp1),tmpancestry_specific,tcga_variant_info[match(rownames(tmp1),tcga_variant_info$HGVSg1),c('HGVSp','Function','ACMG','Symbol_Merge')],stringsAsFactors = F)
    colnames(tmp2)[1:6]<-c('TCGA_AC','Gnomad_AC','TCGA_AF','Gnomad_AF','HGVSg','Specificity')
    tmp2$Label<-paste(tmp2$Symbol_Merge,sapply(strsplit(tmp2$HGVSp,':'),function(x){x[2]}),sep = ':')
    
    tmpfun<-function(x,n){if(length(x)>n){return(which(x>=sort(x,decreasing = T)[n]))}else{return(1:length(x))}}#for function
    if(sum(tmp2$Specificity=='Yes')>0){
      tmptmp2<-tmp2[tmp2$Specificity=='Yes',]
      tmptmp3<-unique(rbind(tmp2[union(tmpfun(tmp2$TCGA_AC,1),tmpfun(tmp2$Gnomad_AC,1)),],tmptmp2[union(tmpfun(tmptmp2$TCGA_AC,1),tmpfun(tmptmp2$Gnomad_AC,1)),]))
    }else{tmptmp3<-unique(tmp2[union(tmpfun(tmp2$TCGA_AC,1),tmpfun(tmp2$Gnomad_AC,1)),])}#for if
    
    m[i,'Num_plot']<-nrow(tmp2); m[i,'Num_label']<-nrow(tmptmp3)
    m[i,'Num_spe']<-sum(tmp2$Specificity=='Yes'); m[i,'Num_not_spe']<-sum(tmp2$Specificity=='No')
    
    if(m[i,'Num_plot']>10){
      tmp<-cor.test(tmp2$TCGA_AF,tmp2$Gnomad_AF,method = 'pearson')
      m[i,'Freq_pearson']<-tmp$estimate; m[i,'Freq_pearson_p']<-tmp$p.value
      tmp<-cor.test(tmp2$TCGA_AF,tmp2$Gnomad_AF,method = 'spearman')
      m[i,'Freq_spearman']<-tmp$estimate; m[i,'Freq_spearman_p']<-tmp$p.value
      tmp<-cor.test(tmp2$TCGA_AC,tmp2$Gnomad_AC,method = 'pearson')
      m[i,'Count_pearson']<-tmp$estimate; m[i,'Count_pearson_p']<-tmp$p.value
      tmp<-cor.test(tmp2$TCGA_AC,tmp2$Gnomad_AC,method = 'spearman')
      m[i,'Count_spearman']<-tmp$estimate; m[i,'Count_spearman_p']<-tmp$p.value
    }#for if    
  }#for if
}#for i
m1<-m[c('Afr','Amr','Eur','Eas','Sas','Oth'),]
write.table(m1,'analysis/data/gnomad/ancestry_variant_distribution/statistic/statistic.txt',row.names = T,col.names = T,quote = F,sep = '\t')

