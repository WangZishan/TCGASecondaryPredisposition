rm(list = ls()); tmp<-lapply(c('ggplot2','reshape2'), library, character.only=T)
setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition')
indir_ase<-outdir<-file.path(getwd(),'analysis/allele_expression/SPV')

### load in related information
# ase analysis info
ase<-read.table(file.path(indir_ase,'ASE_basic-allele_count.txt'),header = T,sep = '\t',stringsAsFactors = F)
ase<-ase[ase$class%in%c('Not ASE','Suggestive','Significant'),]
# ancestry info
ancestry<-read.table('analysis/data/ancestry_info/ancestry_patient.txt',header = T,sep = '\t',stringsAsFactors = F)
# clinical info
clin_complete = read.table('../Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/clinical_PANCAN_patient_with_followup.tsv', header=T, quote = '', sep='\t', fill =T, stringsAsFactors=FALSE)

### merge ancestry information and cancer type info
ase$ancestry<-ancestry[match(ase$bcr_patient_barcode,ancestry$bcr_patient_barcode),'ancestry']
ase$Cancer_Type<-clin_complete[match(ase$bcr_patient_barcode,clin_complete[,'bcr_patient_barcode']),'acronym']
ase$Chromosome<-paste('Chr',ase$Chromosome,sep = '')

### SPVs count/ratio across different levels of a specific aspect of interest
###------------------------------------------------------------------------------------------------------------------------------
### change the aspect name to switch between different aspects of interest
###------------------------------------------------------------------------------------------------------------------------------
aspect<-'Function'
# count
long<-as.data.frame(table(ase[,c('class',aspect)]),stringsAsFactors = F); colnames(long)<-c('class','aspect','Freq')
long$class<-factor(long$class,levels = c('Insufficient read count data','Not ASE','Suggestive','Significant'))
if(aspect=='ancestry'){long$aspect<-factor(long$aspect,levels = c('Amr','Afr','Eur','Eas','Sas','Mix','Oth'))}
ASE_status_color<-c('Significant'=rgb(97,89,164,maxColorValue = 255),'Suggestive'=rgb(206,74,8,maxColorValue = 255),'Not ASE'=rgb(29,143,100,maxColorValue = 255),'Insufficient read count data'='gray70')
p<-ggplot(data = long,aes(x=aspect,y=Freq,fill=class))+geom_bar(stat = 'identity')
p<-p+scale_fill_manual(values = ASE_status_color)
p<-p+theme_bw()+theme(axis.text.x = element_text(size = 18,angle = 90,hjust = 1,vjust = 0.4),axis.text.y = element_text(size = 18),
                      axis.title = element_text(size = 20),legend.text = element_text(size = 18),legend.title = element_text(size = 20))
p
ggsave(file.path(outdir,paste('ASE_MultiVars-',aspect,'_Count.pdf',sep = '')),h=4,w=9,useDingbat=F)
# ratio
count_wide<-dcast(long,value.var = 'Freq',formula = class~aspect)
write.table(count_wide,file.path(outdir,paste('ASE_MultiVars-',aspect,'_Count.txt',sep = '')),row.names = F,col.names = T,quote = F,sep = '\t')

tmp<-apply(count_wide[,2:ncol(count_wide)],2,sum); long$ratio<-long$Freq/tmp[as.character(long$aspect)]
p<-ggplot(data = long,aes(x=aspect,y=ratio,fill=class))+geom_bar(stat = 'identity')
p<-p+scale_fill_manual(values = ASE_status_color)
p<-p+theme_bw()+theme(axis.text.x = element_text(size = 18,angle = 90,hjust = 1,vjust = 0.4),axis.text.y = element_text(size = 18),
                      axis.title = element_text(size = 20),legend.text = element_text(size = 18),legend.title = element_text(size = 20))
p
ggsave(file.path(outdir,paste('ASE_MultiVars-',aspect,'_Proportion.pdf',sep = '')),h=4,w=9,useDingbat=F)
ratio_wide<-dcast(long,value.var = 'ratio',formula = class~aspect)
write.table(ratio_wide,file.path(outdir,paste('ASE_MultiVars-',aspect,'_Proportion.txt',sep = '')),row.names = F,col.names = T,quote = F,sep = '\t')
ratio_wide1<-data.frame(count_wide$class,count_wide[,2:ncol(count_wide)]/apply(count_wide[,2:ncol(count_wide)],1,sum)); colnames(ratio_wide1)[1]<-'class'
write.table(ratio_wide1,file = file.path(outdir,paste('ASE_MultiVars-',aspect,'_Proportion1.txt',sep = '')),row.names = F,col.names = T,quote = F,sep = '\t')

### Permutation analysis
sig<-apply(ratio_wide[match('Significant',ratio_wide$class),2:ncol(ratio_wide)],2,sum)
sig_sug<-apply(ratio_wide[match(c('Suggestive','Significant'),ratio_wide$class),2:ncol(ratio_wide)],2,sum)
rand<-10000; vars<-sort(unique(ase[,aspect]))
m<-matrix(0,nrow = length(vars),ncol = 2,dimnames = list(vars,c('pval_sig','pval_sig_sug')))
for (i in 1:rand) {
  tmpclass<-ase$class; tmptype<-as.character(sample(ase[,aspect],nrow(ase)))
  tmpase<-data.frame(tmpclass,tmptype,stringsAsFactors = F)
  tmpsta<-dcast(as.data.frame(table(tmpase),stringsAsFactors = F),value.var = 'Freq',formula =  tmptype ~ tmpclass)
  
  tmptype<-as.character(tmpsta$tmptype);  tmpsta<-data.frame(tmptype,tmpsta[,-1]/apply(tmpsta[,-1],1,sum))
  tmpsig<-tmpsta$Significant;  tmpsig_sug<-apply(tmpsta[,c('Significant','Suggestive')], 1, sum)
  names(tmpsig)<-names(tmpsig_sug)<-tmptype
  
  for (j in 1:length(tmptype)) {
    tmptmptype<-tmptype[j]
    if(tmpsig[tmptmptype]>sig[tmptmptype]){m[tmptmptype,'pval_sig']<-m[tmptmptype,'pval_sig']+1}
    if(tmpsig_sug[tmptmptype]>sig_sug[tmptmptype]){m[tmptmptype,'pval_sig_sug']<-m[tmptmptype,'pval_sig_sug']+1}
  }#for j
}#for i
m<-m/rand; m
fdr_sig<-p.adjust(m[,'pval_sig'],method = 'BH'); fdr_sig_sug<-p.adjust(m[,'pval_sig_sug'],method = 'BH')
m<-data.frame(m,fdr_sig,fdr_sig_sug,stringsAsFactors = F)
write.table(m,file.path(outdir,paste('ASE_MultiVars-',aspect,'_Significance.txt',sep = '')),row.names = T,col.names = T,quote = F,sep = '\t')

