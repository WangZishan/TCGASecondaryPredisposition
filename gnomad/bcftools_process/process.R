rm(list = ls()); setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition')
firstup<-function(x){substr(x,1,1)<-toupper(substr(x,1,1));return(x)}
outdir<-'C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition/analysis/data/gnomad/bcftools_process/process'

### Load in related information
# variant information
variant<-read.table('analysis/data/variant_data/variants_pca_frq_nonCancer_pathogenic.txt',sep = '\t',stringsAsFactors = F,header = T,quote = '')
# gnomad variant information
gnomad_acmg_equal<-read.table('analysis/data/gnomad/bcftools/VariantEqualPos.vcf',sep = '\t',stringsAsFactors = F)
gnomad_acmg_notequal<-read.table('analysis/data/gnomad/bcftools/VariantNotequalPos.vcf',sep = '\t',stringsAsFactors = F)
gnomad_acmg <- unique(rbind(gnomad_acmg_equal,gnomad_acmg_notequal))
colnames(gnomad_acmg)<-c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO')
gnomad_acmg<-gnomad_acmg[gnomad_acmg$FILTER=='PASS',]
gnomad_acmg$HGVSg1 <- NA
for (i in 1:nrow(gnomad_acmg)) {
  tmp1 <- nchar(gnomad_acmg$REF[i]); tmp2 <- nchar(gnomad_acmg$ALT[i])
  if(tmp1==1 & tmp2==1){
    gnomad_acmg$HGVSg1[i] <- paste(gnomad_acmg$CHROM[i],':g.',gnomad_acmg$POS[i],gnomad_acmg$REF[i],'>',gnomad_acmg$ALT[i],sep = '')
  }else if(tmp1>1 & tmp2==1){
    gnomad_acmg$HGVSg1[i] <- paste(gnomad_acmg$CHROM[i],':g.',gnomad_acmg$POS[i]+1,'_',gnomad_acmg$POS[i]+nchar(gnomad_acmg$REF[i])-1,'del',substr(gnomad_acmg$REF[i],2,nchar(gnomad_acmg$REF[i])),sep = '')
  }else if(tmp1==1 & tmp2>1){
    gnomad_acmg$HGVSg1[i] <- paste(gnomad_acmg$CHROM[i],':g.',gnomad_acmg$POS[i],'_',gnomad_acmg$POS[i]+1,'ins',substr(gnomad_acmg$ALT[i],2,nchar(gnomad_acmg$ALT[i])),sep = '')
  }
}# for
gnomad_acmg <- gnomad_acmg[,c('HGVSg1','INFO')]

### Extract interested information for our variants from gnomad
gnomad_acmg1 <- unique(merge(gnomad_acmg,unique(variant[,c('HGVSg1','Chromosome','Start','Stop','Reference','ALT_clean','Function','ACMG','Symbol_Merge')])))
ancestries<-c('amr','afr','fin','nfe','asj','eas','sas','oth')
vars<-c('AC','AN','non_cancer_AC','non_cancer_AN')
for (i in 1:length(vars)) {
  minfo<-as.data.frame(matrix(nrow = nrow(gnomad_acmg1),ncol = length(ancestries),dimnames = list(NULL,paste(vars[i],ancestries,sep = '_'))),stringsAsFactors = F)
  for (j in 1:nrow(gnomad_acmg1)) {
      tmp<-unlist(strsplit(gnomad_acmg1[j,'INFO'],';'))
      for (k in 1:length(ancestries)) {
        tmptmp<-tmp[setdiff(grep(paste(vars[i],'_',ancestries[k],'=',sep = ''),tmp),grep(paste('_',vars[i],'_',ancestries[k],'=',sep = ''),tmp))]
        minfo[j,k]<-unlist(strsplit(tmptmp,'='))[2]
      }#for k
  }#for j
  minfo<-apply(minfo,2,as.numeric)
  minfo_asian<-apply(minfo[,paste(vars[i],c('eas','sas'),sep = '_')],1,function(x){if(any(is.na(x))){return(NA)}else{sum(x)}})
  minfo_eur<-apply(minfo[,paste(vars[i],c('fin','nfe'),sep = '_')],1,function(x){if(any(is.na(x))){return(NA)}else{sum(x)}})
  minfo1<-data.frame(minfo,minfo_asian,minfo_eur);  colnames(minfo1)<-paste(firstup(c(ancestries,'asian','eur')),vars[i],sep = '_')
  minfo2<-minfo1[,paste(firstup(c('amr','afr','fin','nfe','asj','eur','eas','sas','asian','oth')),vars[i],sep = '_')]
  assign(vars[i],minfo2)
}#for i
AF<-AC/AN; colnames(AF)<-gsub('AC','AF',colnames(AF))
non_cancer_AF<-non_cancer_AC/non_cancer_AN; colnames(non_cancer_AF)<-gsub('AC','AF',colnames(non_cancer_AF))
gnomad_acmg2<-data.frame(gnomad_acmg1[,setdiff(colnames(gnomad_acmg1),'INFO')],AC,AN,AF,non_cancer_AC,non_cancer_AN,non_cancer_AF)
write.table(gnomad_acmg2,file.path(outdir,'gnomad_variant_info.txt'),row.names = F,col.names = T,quote = F,sep = '\t')
vars_interest<-c('AC','AN','AF','non_cancer_AC','non_cancer_AN','non_cancer_AF')
for (i in vars_interest) {write.table(get(i),file.path(outdir,paste(i,'txt',sep = '.')),row.names = F,col.names = T,quote = F,sep = '\t')}#for i

#collect interested information for each ACMG status (pathogenic and likely pathogenic)
acmg_status<-unique(gnomad_acmg2$ACMG)
for (i in 1:length(vars_interest)) {
  l<-list()
  for (j in 1:length(acmg_status)) {
    tmptmp<-get(vars_interest[i])[gnomad_acmg2$ACMG==acmg_status[j],]
    if(length(grep('AN',vars_interest[i]))==1){l[[acmg_status[j]]]<-colMeans(tmptmp)}else{l[[acmg_status[j]]]<-colSums(tmptmp)}
  }#for j
  assign(paste(vars_interest[i],'acmg_status',sep = '_'),do.call(rbind,l))
}#for i
AF1_acmg_status<-AC_acmg_status/AN_acmg_status; colnames(AF1_acmg_status)<-gsub('AC','AF1',colnames(AF1_acmg_status))
non_cancer_AF1_acmg_status<-non_cancer_AC_acmg_status/non_cancer_AN_acmg_status; colnames(non_cancer_AF1_acmg_status)<-gsub('AC','AF1',colnames(non_cancer_AF1_acmg_status))
for (i in c(vars_interest,'AF1','non_cancer_AF1')) {write.table(get(paste(i,'acmg_status',sep = '_')),file.path(outdir,paste(i,'acmg_status.txt',sep = '_')),row.names = T,col.names = T,quote = F,sep = '\t')}#for i

#collect interested information for each gene
genes<-unique(gnomad_acmg2$Symbol_Merge)
for (i in 1:length(vars_interest)) {
  l<-list()
  for (j in 1:length(genes)) {
    tmptmp<-get(vars_interest[i])[gnomad_acmg2$Symbol_Merge==genes[j],]
    if(length(grep('AN',vars_interest[i]))==1){l[[genes[j]]]<-colMeans(tmptmp)}else{l[[genes[j]]]<-colSums(tmptmp)}
  }#for j
  assign(paste(vars_interest[i],'gene',sep = '_'),do.call(rbind,l))
}#for i
AF1_gene<-AC_gene/AN_gene; colnames(AF1_gene)<-gsub('AC','AF1',colnames(AF1_gene))
non_cancer_AF1_gene<-non_cancer_AC_gene/non_cancer_AN_gene; colnames(non_cancer_AF1_gene)<-gsub('AC','AF1',colnames(non_cancer_AF1_gene))
for (i in c(vars_interest,'AF1','non_cancer_AF1')) {write.table(get(paste(i,'gene',sep = '_')),file.path(outdir,paste(i,'gene.txt',sep = '_')),row.names = T,col.names = T,quote = F,sep = '\t')}#for i

