###--------------------------------------------------------------------------------------------------------------------------------------------------
### The difference between sig and sig sug is to change filter condition and output filename
###--------------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls()); tmp<-lapply(c('ggplot2','ggbeeswarm'), library, character.only=T)
setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition')
outdir<-'C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition/analysis/VariantImpactOnExp/exp_percentile/'
ASE_status_color<-c('Significant'=rgb(97,89,164,maxColorValue = 255),'Suggestive'=rgb(206,74,8,maxColorValue = 255),'Not ASE'=rgb(29,143,100,maxColorValue = 255),'Insufficient read count data'='gray70')
Variant_type_color<-c('frameshift'='#66C2A5','missense'='#FC8D62','splice-variant'='#8DA0CB','splice-acceptor'='#A6D854','splice-donor'='#33A02C',
                      'start-lost'='#E78AC3','stop-gain'='#A6761D','stop-lost'='#E5C494','synonymous'='#FFD92F')
mydiscrete<-function(x,binwidth=0.01,minvalue=0,maxvalue=1){
  val<-c(seq(minvalue,maxvalue,by = binwidth),1)
  x1<-x
  for (i in 1:length(x)) {
    tmpx<-x[i]
    tmpx_dis<-abs(tmpx-val)
    min_ind<-which(tmpx_dis==min(tmpx_dis))[1]
    x1[i]<-val[min_ind]
  }#for x
  return(x1)
}#for function

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
genes_enriched_ase<-tmpgenes_enriched_ase$Symbol_Merge[tmpgenes_enriched_ase$Classification%in%c('Significant','Suggestive')]
#Genes whose expression impacted by variant using linear regression model
files<-c('individual_cancer','pan_cancer')
for (i in 1:length(files)) {
  tmp<-read.table(paste('analysis/VariantImpactOnExp/exp_percentile/VariantImpactOnExp-',files[i],'.txt',sep = ''),header = T,sep = '\t',stringsAsFactors = F)
  if(length(grep('individual',files[i]))==1){
    tmp1<-tmp[!is.na(tmp$P_Chi1)&(tmp$P_Chi1<0.15),]; tmp2<-tmp1[order(tmp1$P_Chi1),]
  }else{tmp1<-tmp[!is.na(tmp$FDR1)&(tmp$FDR1<0.15),]; tmp2<-tmp1[order(tmp1$FDR1),]}
  assign(files[i],tmp2)
}#for i

###Combine ase and variant type information to plot for percentile expression of genes, which is either enriched with ase or impacted by variant
#define a function to plot this kind of data structure
plot_percentile<-function(genes,genes_variant,ase_info,express,dir_plot,fn){
  tmpgenes_variant<-genes_variant[genes_variant$Symbol_Merge%in%genes,c('Symbol_Merge','HGVSg1','bcr_patient_barcode','Sample','Function')]
  tmpase_info<-ase_info[,c('Symbol_Merge','HGVSg1','bcr_patient_barcode','Sample','class')]
  tmpmerge1<-merge(tmpgenes_variant,tmpase_info,all.x = T); tmpmerge1$class[is.na(tmpmerge1$class)]<-'Insufficient read count data'
  expression<-apply(tmpmerge1[,c('Symbol_Merge','bcr_patient_barcode')],1,function(x){express[x[1],x[2]]})
  tmpmerge<-data.frame(tmpmerge1,expression,stringsAsFactors = F)
  tmpmerge$Symbol_Merge<-factor(tmpmerge$Symbol_Merge,levels=genes)
  
  p<-ggplot(tmpmerge,aes(x=Symbol_Merge,y=expression,colour=class,fill=Function,stroke=3.5))
  p<-p+geom_dotplot(dotsize=2,binwidth=.01, binaxis= 'y',stackdir ='centerwhole')
  p<-p+scale_color_manual(name='ASE Status',values = ASE_status_color)
  p<-p+scale_fill_manual(name='Variant Type',values = Variant_type_color)
  p<-p+labs(x='Genes',y='Expression')
  p<-p+theme_minimal()+theme(axis.title = element_text(size = 22),axis.text.y = element_text(size = 20),
             axis.text.x = element_text(size = 20,angle = 90,vjust = 0.5,hjust = 1),
             legend.title = element_text(size = 22),legend.text = element_text(size = 20),
             legend.key.size = unit(1.3,'cm'),legend.box.background = element_blank(),legend.position = 'bottom',
             panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
  p
  ggsave(h=9,w=2+length(genes)/3,useDingbat=F,filename = file.path(dir_plot,fn))
}#for function
#Plot for the percentile expression of genes, which is enriched with ase
plot_percentile(genes = genes_enriched_ase,
                genes_variant = variant_cancer,
                ase_info = ase_info,
                express = exp_cancer,
                dir_plot = outdir,
                fn = 'SugSig_PlotPercentileExp-Combine-Genes_enriched_with_ASE-Cancer_Exp.pdf')
#Plot for the percentile expression of genes impacted by variant using linear regression model
for (i in 1:length(files)) {
  tmp<-get(files[i])
  if(nrow(tmp)!=0){
    plot_percentile(genes = unique(tmp$gene),
                    genes_variant = variant_cancer,
                    ase_info = ase_info,
                    express = exp_cancer,
                    dir_plot = outdir,
                    fn = paste('SugSig_PlotPercentileExp-','Combine-',files[i],'.pdf',sep = ''))
  }#for if
}#for i
#Plot for the percentile expression of genes in the union set of enriched ase and impacted by linear regression model
for (i in 1:length(files)) {
  tmp<-get(files[i])
  if(nrow(tmp)!=0){
    tmpgene_inter<-sort(intersect(tmp$gene,genes_enriched_ase))
    cat(files[i]); cat('\n'); cat(tmpgene_inter)
    tmpaes_spe<-setdiff(genes_enriched_ase,tmpgene_inter)
    tmpimpact_spe<-setdiff(tmp$gene,tmpgene_inter)
    tmpgene<-c(tmpimpact_spe,tmpgene_inter,tmpaes_spe)
    plot_percentile(genes = tmpgene,
                    genes_variant = variant_cancer,
                    ase_info = ase_info,
                    express = exp_cancer,
                    dir_plot = outdir,
                    fn = paste('SugSig_PlotPercentileExp-','Combine-',files[i],'_union.pdf',sep = ''))
  }#for if
}#for i

###Only use ase information to plot for percentile expression of genes, which is either enriched with ase or impacted by variant
#define a function to plot this kind of data structure
plot_percentile<-function(genes,genes_variant,ase_info,express,dir_plot,fn){
  tmpgenes_variant<-genes_variant[genes_variant$Symbol_Merge%in%genes,c('Symbol_Merge','HGVSg1','bcr_patient_barcode','Sample','Function')]
  tmpase_info<-ase_info[,c('Symbol_Merge','HGVSg1','bcr_patient_barcode','Sample','class')]
  tmpmerge1<-merge(tmpgenes_variant,tmpase_info,all.x = T);tmpmerge1$class[is.na(tmpmerge1$class)]<-'Insufficient read count data'
  expression<-apply(tmpmerge1[,c('Symbol_Merge','bcr_patient_barcode')],1,function(x){express[x[1],x[2]]})
  tmpmerge<-data.frame(tmpmerge1,expression,stringsAsFactors = F)
  tmpmerge$Symbol_Merge<-factor(tmpmerge$Symbol_Merge,levels=genes)
  tmpmerge$expression1<-mydiscrete(x=tmpmerge$expression,binwidth = 0.01)
  
  p<-ggplot(tmpmerge,aes(x=Symbol_Merge,y=expression1,color=class,stroke=0))
  p<-p+geom_beeswarm(size=2.5)
  p<-p+scale_color_manual(name='ASE Status',values = ASE_status_color)
  p<-p+labs(x='Genes',y='Expression')
  p<-p+theme_minimal()
  p<-p+theme(axis.title = element_text(size = 22),axis.text.y = element_text(size = 20),
             axis.text.x = element_text(size = 20,angle = 90,vjust = 0.5,hjust = 1),
             legend.title = element_text(size = 22),legend.text = element_text(size = 20),
             legend.key.size = unit(1.3,'cm'),legend.box.background = element_blank(),legend.position = 'bottom',
             panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
  p
  ggsave(h=9,w=2+length(genes)/3,useDingbat=F,filename = file.path(dir_plot,fn))
}#for function
#Plot for percentile expression of genes, which is enriched with ase
plot_percentile(genes = genes_enriched_ase,
                genes_variant = variant_cancer,
                ase_info = ase_info,
                express = exp_cancer,
                dir_plot = outdir,
                fn = 'SugSig_PlotPercentileExp-Only_ase-Genes_enriched_with_ASE-Cancer_Exp.pdf')
#Plot for percentile expression of genes impacted by variant using linear regression model
for (i in 1:length(files)) {
  tmp<-get(files[i])
  if(nrow(tmp)!=0){
    plot_percentile(genes = unique(tmp$gene),
                    genes_variant = variant_cancer,
                    ase_info = ase_info,
                    express = exp_cancer,
                    dir_plot = outdir,
                    fn = paste('SugSig_PlotPercentileExp-','Only_ase-',files[i],'.pdf',sep = ''))
  }#for if
}#for i
#Plot for percentile expression of genes in the union set of enriched ase and impacted by linear regression model
for (i in 1:length(files)) {
  tmp<-get(files[i])
  if(nrow(tmp)!=0){
    tmpgene_inter<-sort(intersect(tmp$gene,genes_enriched_ase))
    cat(files[i])
    cat('\n')
    cat(tmpgene_inter)
    tmpaes_spe<-setdiff(genes_enriched_ase,tmpgene_inter)
    tmpimpact_spe<-setdiff(tmp$gene,tmpgene_inter)
    tmpgene<-c(tmpimpact_spe,tmpgene_inter,tmpaes_spe)
    plot_percentile(genes = tmpgene,
                    genes_variant = variant_cancer,
                    ase_info = ase_info,
                    express = exp_cancer,
                    dir_plot = outdir,
                    fn = paste('SugSig_PlotPercentileExp-','Only_ase-',files[i],'_union.pdf',sep = ''))
  }#for if
}#for i

###only use variant type information to plot for percentile expression of genes, which is either enriched with ase or impacted by variant
#define a function to plot this kind of data structure
plot_percentile<-function(genes,genes_variant,ase_info,express,dir_plot,fn){
  tmpgenes_variant<-genes_variant[genes_variant$Symbol_Merge%in%genes,c('Symbol_Merge','HGVSg1','bcr_patient_barcode','Sample','Function')]
  tmpase_info<-ase_info[,c('Symbol_Merge','HGVSg1','bcr_patient_barcode','Sample','class')]
  tmpmerge1<-merge(tmpgenes_variant,tmpase_info,all.x = T);tmpmerge1$class[is.na(tmpmerge1$class)]<-'Insufficient read count data'
  expression<-apply(tmpmerge1[,c('Symbol_Merge','bcr_patient_barcode')],1,function(x){express[x[1],x[2]]})
  tmpmerge<-data.frame(tmpmerge1,expression,stringsAsFactors = F)
  tmpmerge$Symbol_Merge<-factor(tmpmerge$Symbol_Merge,levels=genes)
  tmpmerge$expression1<-mydiscrete(x=tmpmerge$expression,binwidth = 0.01)
  
  p<-ggplot(tmpmerge,aes(x=Symbol_Merge,y=expression1,color=Function,stroke=0))
  p<-p+geom_beeswarm(size=2.5)
  p<-p+scale_color_manual(name='Variant Type',values = Variant_type_color)
  p<-p+labs(x='Genes',y='Expression')
  p<-p+theme_minimal()
  p<-p+theme(axis.title = element_text(size = 22),axis.text.y = element_text(size = 20),
             axis.text.x = element_text(size = 20,angle = 90,vjust = 0.5,hjust = 1),
             legend.title = element_text(size = 22),legend.text = element_text(size = 20),
             legend.key.size = unit(1.3,'cm'),legend.box.background = element_blank(),legend.position = 'bottom',
             panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
  p
  ggsave(h=9,w=2+length(genes)/3,useDingbat=F,filename = file.path(dir_plot,fn))
}#for function
#Plot for percentile expression of genes, which is enriched with ase
plot_percentile(genes = genes_enriched_ase,
                genes_variant = variant_cancer,
                ase_info = ase_info,
                express = exp_cancer,
                dir_plot = outdir,
                fn = 'SugSig_PlotPercentileExp-Only_variant_type-Genes_enriched_with_ASE-Cancer_Exp.pdf')
#Plot for percentile expression of genes impacted by variant using linear regression model
for (i in 1:length(files)) {
  tmp<-get(files[i])
  if(nrow(tmp)!=0){
    plot_percentile(genes = unique(tmp$gene),
                    genes_variant = variant_cancer,
                    ase_info = ase_info,
                    express = exp_cancer,
                    dir_plot = outdir,
                    fn = paste('SugSig_PlotPercentileExp-','Only_variant_type-',files[i],'.pdf',sep = ''))
  }#for if
}#for i
#Plot for percentile expression of genes in the union set of enriched ase and impacted by linear regression model
for (i in 1:length(files)) {
  tmp<-get(files[i])
  if(nrow(tmp)!=0){
    tmpgene_inter<-sort(intersect(tmp$gene,genes_enriched_ase))
    cat(files[i])
    cat('\n')
    cat(tmpgene_inter)
    tmpaes_spe<-setdiff(genes_enriched_ase,tmpgene_inter)
    tmpimpact_spe<-setdiff(tmp$gene,tmpgene_inter)
    tmpgene<-c(tmpimpact_spe,tmpgene_inter,tmpaes_spe)
    plot_percentile(genes = tmpgene,
                    genes_variant = variant_cancer,
                    ase_info = ase_info,
                    express = exp_cancer,
                    dir_plot = outdir,
                    fn = paste('SugSig_PlotPercentileExp-','Only_variant_type-',files[i],'_union.pdf',sep = ''))
  }#for if
}#for i

