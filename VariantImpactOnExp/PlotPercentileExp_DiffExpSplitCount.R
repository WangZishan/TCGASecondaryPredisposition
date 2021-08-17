###--------------------------------------------------------------------------------------------------------------------------------------------------
### The difference between sig and sig sug is to change filter condition and output filename
###--------------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls()); setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition')
sta_bin<-function(x){x<-x[!is.na(x)];return(c(sum(x<0.25),sum(x<0.5&x>=0.25),sum(x<0.75&x>=0.5),sum(x>=0.75)))}#for function

### Load in related data
# expression data of variant genes
exp_normal<-read.table('analysis/VariantImpactOnExp/data/expression_of_variant_genes_preprocess/normal_expression_of_variant_genes-pc-clinical-ecdf.txt',row.names = 1,header = T,sep = '\t',stringsAsFactors = F)
exp_cancer<-read.table('analysis/VariantImpactOnExp/data/expression_of_variant_genes_preprocess/cancer_expression_of_variant_genes-pc-clinical-ecdf.txt',row.names = 1,header = T,sep = '\t',stringsAsFactors = F)
colnames(exp_normal)<-gsub('\\.','-',colnames(exp_normal)); colnames(exp_cancer)<-gsub('\\.','-',colnames(exp_cancer))
exp_normal_sample<-read.table('analysis/VariantImpactOnExp/data/expression_of_variant_genes_preprocess/normal_expression_of_variant_genes-pc-clinical-sampleinfo.txt',header = T,sep = '\t',stringsAsFactors = F)
exp_cancer_sample<-read.table('analysis/VariantImpactOnExp/data/expression_of_variant_genes_preprocess/cancer_expression_of_variant_genes-pc-clinical-sampleinfo.txt',header = T,sep = '\t',stringsAsFactors = F)
#variant data
variant<-read.table('analysis/data/variant_data/variants_pca_frq_nonCancer_pathogenic.txt',header = T,sep = '\t',stringsAsFactors = F,quote = '')
variant_cancer<-variant[variant$bcr_patient_barcode%in%colnames(exp_cancer),]
#Genes enriched with ASE
ase_info<-read.table('analysis/allele_expression/SPV/ASE_basic-allele_count.txt',header = T,sep = '\t',stringsAsFactors = F)
tmpgenes_enriched_ase <- read.table('analysis/allele_expression/SPV/ASE_Gene_Enrichment.txt',header = T,sep = '\t',stringsAsFactors = F)
genes_enriched_ase<-tmpgenes_enriched_ase$Symbol_Merge[tmpgenes_enriched_ase$Classification%in%c('Significant')]
#Genes whose expression impacted by variant using linear regression model
files<-c('individual_cancer','pan_cancer')
for (i in 1:length(files)) {
  tmp<-read.table(paste('analysis/VariantImpactOnExp/exp_percentile/VariantImpactOnExp-',files[i],'.txt',sep = ''),header = T,sep = '\t',stringsAsFactors = F)
  if(length(grep('individual',files[i]))==1){
    tmp1<-tmp[!is.na(tmp$P_Chi1)&(tmp$P_Chi1<0.05),]; tmp2<-tmp1[order(tmp1$P_Chi1),]
  }else{tmp1<-tmp[!is.na(tmp$FDR1)&(tmp$FDR1<0.05),]; tmp2<-tmp1[order(tmp1$FDR1),]}
  assign(files[i],tmp2)
}#for i

### Count the number of sample-variants of different expression splits across different ASE status level/variant types for a specific genes set
Gene_VariantImpactOnExp <- pan_cancer$gene; Gene_EnrichedASE <- genes_enriched_ase
Gene_VariantImpactOnExp_EnrichedASE <- union(Gene_VariantImpactOnExp,Gene_EnrichedASE); Gene_All <- unique(variant_cancer$Symbol_Merge)
geneset <- c('Gene_VariantImpactOnExp','Gene_EnrichedASE','Gene_VariantImpactOnExp_EnrichedASE','Gene_All')
for (j in 1:length(geneset)) {
  tmpvariant_cancer<-variant_cancer[variant_cancer$Symbol_Merge%in%get(geneset[j]),c('Symbol_Merge','HGVSg1','bcr_patient_barcode','Sample','Function')]
  tmpase_info<-ase_info[,c('Symbol_Merge','HGVSg1','bcr_patient_barcode','Sample','class')]
  tmpmerge1<-merge(tmpvariant_cancer,tmpase_info,all.x = T); tmpmerge1$class[is.na(tmpmerge1$class)]<-'Insufficient read count data'
  expression<-apply(tmpmerge1[,c('Symbol_Merge','bcr_patient_barcode')],1,function(x){exp_cancer[x[1],x[2]]})
  tmpmerge<-data.frame(tmpmerge1,expression,stringsAsFactors = F)
  
  Total<-sta_bin(tmpmerge$expression); Total<-c(Total,NA,Total/sum(Total))
  ASE_status<-c('Significant','Suggestive','Not ASE','Insufficient read count data')
  m1<-matrix(nrow = length(ASE_status),ncol =4,dimnames = list(ASE_status,c('0-0.25','0.25-0.5','0.5-0.75','0.75-1')))
  for (i in 1:length(ASE_status)) {m1[i,]<-sta_bin(tmpmerge[tmpmerge$class==ASE_status[i],'expression'])}#for i
  m1<-cbind(m1,NA,m1/apply(m1, 1, sum))
  Variant_type<-sort(unique(variant$Function))
  m2<-matrix(nrow = length(Variant_type),ncol =4,dimnames = list(Variant_type,c('0-0.25','0.25-0.5','0.5-0.75','0.75-1')))
  for (i in 1:length(Variant_type)) {m2[i,]<-sta_bin(tmpmerge[tmpmerge$Function==Variant_type[i],'expression'])}#for i
  m2<-cbind(m2,NA,m2/apply(m2, 1, sum))
  m<-rbind(m1,NA,m2,NA,Total);m<-data.frame(rownames(m),m); colnames(m)<-c('Class',colnames(m1))
  writexl::write_xlsx(m,file.path('analysis/VariantImpactOnExp/exp_percentile',paste('Sig_PlotPercentileExp_DiffExpSplitCount-',geneset[j],'.xlsx',sep = '')),format_headers = F)
}# for j


