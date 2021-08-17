###--------------------------------------------------------------------------------------------------------------------------------------------------
### The difference between sig and sig sug is to change filter condition and output filename
###--------------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls()); setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition')

### Load in related information
# expression data of variant genes
exp_normal<-read.table('analysis/VariantImpactOnExp/data/expression_of_variant_genes_preprocess/normal_expression_of_variant_genes-pc-clinical-ecdf.txt',row.names = 1,header = T,sep = '\t',stringsAsFactors = F)
exp_cancer<-read.table('analysis/VariantImpactOnExp/data/expression_of_variant_genes_preprocess/cancer_expression_of_variant_genes-pc-clinical-ecdf.txt',row.names = 1,header = T,sep = '\t',stringsAsFactors = F)
colnames(exp_normal)<-gsub('\\.','-',colnames(exp_normal)); colnames(exp_cancer)<-gsub('\\.','-',colnames(exp_cancer))
exp_normal_sample<-read.table('analysis/VariantImpactOnExp/data/expression_of_variant_genes_preprocess/normal_expression_of_variant_genes-pc-clinical-sampleinfo.txt',header = T,sep = '\t',stringsAsFactors = F)
exp_cancer_sample<-read.table('analysis/VariantImpactOnExp/data/expression_of_variant_genes_preprocess/cancer_expression_of_variant_genes-pc-clinical-sampleinfo.txt',header = T,sep = '\t',stringsAsFactors = F)
#variant data
variant<-read.table('analysis/data/variant_data/variants_pca_frq_nonCancer_pathogenic.txt',header = T,sep = '\t',stringsAsFactors = F,quote = '')
variant_cancer<-variant[variant$bcr_patient_barcode%in%colnames(exp_cancer),]; variant_normal<-variant[variant$bcr_patient_barcode%in%colnames(exp_normal),]
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

### Extracted information of interested Gene Set
genes_inter<-sort(intersect(pan_cancer$gene,genes_enriched_ase))
genes<-c(setdiff(pan_cancer$gene,genes_inter),genes_inter,setdiff(genes_enriched_ase,genes_inter))
info<-variant[variant$Symbol_Merge%in%genes & variant$bcr_patient_barcode%in%colnames(exp_cancer),c('HGVSg1','Start','Stop','Reference','Alternate','ALT_clean','bcr_patient_barcode','Symbol_Merge','HGVSp','Function')]

ase_info1<-ase_info[,c('HGVSg1','bcr_patient_barcode','class')]
info1<-merge(info,ase_info1,all.x = T)
for (i in 1:nrow(info1)) {info1[i,'ecdf']<-exp_cancer[as.character(info1[i,'Symbol_Merge']),as.character(info1[i,'bcr_patient_barcode'])]}#for i
info1$Symbol_Merge<-factor(info1$Symbol_Merge,levels = genes)
info1<-info1[order(info1$Symbol_Merge,info1$ecdf),]
write.table(info1,'analysis/VariantImpactOnExp/exp_percentile/Sig_PlotPercentileExp_DiffExpSplitCount_GeneInfo.txt',row.names = F,col.names = T,quote = F,sep = '\t')
