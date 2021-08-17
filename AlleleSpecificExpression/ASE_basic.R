###--------------------------------------------------------------------------------------------------------------------------------------------------
### Also used to calculate ASE of more SPVs, which needs to change input file name, output dir and the filter of nonconsitency btw ALT_clean and Alternate
###--------------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls()); tmp<-lapply(c('ggplot2','scales','statmod'), library, character.only=T)
setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition'); source('analysis/functions/general_function.R')
outdir<-file.path(getwd(),'analysis/allele_expression/SPV/')

### load in related info
# variant info
variant<-read.table('analysis/data/variant_data/variants_pca_frq_nonCancer_pathogenic.txt',header = T,sep = '\t',quote = '',stringsAsFactors = F)
# allell count info
allele_count<-read.table('data/variants_pca_frq_rare_all_bamreadcounts.txt',header = T,sep = '\t',stringsAsFactors = F,quote = '')
allele_count<-unique(allele_count[,c('HGVSg','Sample','bcr_patient_barcode','Reference_readcounts','Alternate_readcounts','X3')])
colnames(allele_count)[match(c('HGVSg','Sample','bcr_patient_barcode'),colnames(allele_count))]<-c('HGVSg1','Sample_Count','bcr_patient_barcode_Count')

### preprocess loaded info
#merge varaint data and allele data; select interested allele readcount data of overlapped sample
allele_count<-merge(allele_count,variant)
allele_count<-allele_count[allele_count$Sample==allele_count$Sample_Count,]
allele_count<-unique(allele_count[,c('HGVSg1','Reference','Alternate','ALT_clean','Chromosome','Start','Stop','Sample','bcr_patient_barcode','Function','Symbol_Merge','ACMG','Reference_readcounts','Alternate_readcounts','X3')])
#normalize the data type of readcount
allele_count$X3[is.na(allele_count$X3)]<-NA;allele_count$X3[allele_count$X3=='-' & (!is.na(allele_count$X3))]<-NA
allele_count$Reference_readcounts[is.na(allele_count$Reference_readcounts)]<-NA;allele_count$Reference_readcounts[allele_count$Reference_readcounts=='-']<-NA
allele_count$Alternate_readcounts[is.na(allele_count$Alternate_readcounts)]<-NA;allele_count$Alternate_readcounts[allele_count$Alternate_readcounts=='-']<-NA
allele_count$X3<-as.integer(allele_count$X3);allele_count$Reference_readcounts<-as.integer(allele_count$Reference_readcounts);allele_count$Alternate_readcounts<-as.integer(allele_count$Alternate_readcounts)
allele_count<-allele_count[!is.na(allele_count$Reference_readcounts) & !is.na(allele_count$Alternate_readcounts),]
#compare sequencing coverage with sum of ref readcount and alternate readcount
tmpa<-allele_count[!is.na(allele_count$X3),]
tmpa<-tmpa[tmpa$X3!=(tmpa$Reference_readcounts+tmpa$Alternate_readcounts),]
table(tmpa$X3>(tmpa$Reference_readcounts+tmpa$Alternate_readcounts))

### ASE analysis of each sample-variant by using binom.test
allele_count1<-unique(allele_count[,setdiff(colnames(allele_count),c('Reference_readcounts','Alternate_readcounts','X3'))])
allele_count1_tmpcolnames<-colnames(allele_count1)
var_info<-c('Alternate_readcounts','Reference_readcounts','X3','pval_content','pval_min',
            'sum_max_num','Alternate_readcounts_sum_max','Reference_readcounts_sum_max','pval_sum_max',
            'x3_max_num','Alternate_readcounts_x3_max','Reference_readcounts_x3_max','pval_x3_max')
for (i in 1:length(var_info)) {allele_count1[,var_info[i]]<-NA}
for (i in 1:nrow(allele_count1)) {
  tmpallele_count<-merge(allele_count1[i,allele_count1_tmpcolnames],allele_count)
  allele_count1[i,'Reference_readcounts']<-paste(tmpallele_count$Reference_readcounts,collapse = ',')
  allele_count1[i,'Alternate_readcounts']<-paste(tmpallele_count$Alternate_readcounts,collapse = ',')
  allele_count1[i,'X3']<-paste(tmpallele_count$X3,collapse = ',')
  
  tmp<-tmpallele_count[,c('Reference_readcounts','Alternate_readcounts')]
  tmppval<-apply(tmp,1,function(x){tmptmppval<-binom.test(c(x[1],x[2]),p=0.5,alternative = 'greater')$p.value;return(tmptmppval)})
  allele_count1[i,'pval_content']<-paste(tmppval,collapse = ',')
  tmppval1<-tmppval[!is.na(tmppval)]; if(length(tmppval1)!=0){allele_count1[i,'pval_min']<-min(tmppval1)}#for if
  
  tmpsum<-apply(tmp, 1, sum); tmp1<-tmp[tmpsum==max(tmpsum),]
  allele_count1[i,'Reference_readcounts_sum_max']<-round2(mean(tmp1[,'Reference_readcounts']),0)
  allele_count1[i,'Alternate_readcounts_sum_max']<-round2(mean(tmp1[,'Alternate_readcounts']),0)
  allele_count1[i,'pval_sum_max']<-binom.test(c(allele_count1[i,'Reference_readcounts_sum_max'],allele_count1[i,'Alternate_readcounts_sum_max']),p=0.5,alternative = 'greater')$p.value
  allele_count1[i,'sum_max_num']<-nrow(tmp1)

  tmpind<-!is.na(tmpallele_count$X3)
  if(sum(tmpind)!=0){
    tmpallele_count_nona<-tmpallele_count[tmpind,]
    tmpallele_count_nona1<-tmpallele_count_nona[tmpallele_count_nona$X3==max(tmpallele_count_nona$X3),]
    allele_count1[i,'Reference_readcounts_x3_max']<-round2(mean(tmpallele_count_nona1[,'Reference_readcounts']),0)
    allele_count1[i,'Alternate_readcounts_x3_max']<-round2(mean(tmpallele_count_nona1[,'Alternate_readcounts']),0)
    allele_count1[i,'pval_x3_max']<-binom.test(c(allele_count1[i,'Reference_readcounts_x3_max'],allele_count1[i,'Alternate_readcounts_x3_max']),p=0.5,alternative = 'greater')$p.value
    allele_count1[i,'x3_max_num']<-nrow(tmpallele_count_nona1)
  }#for if
}#for i

### process ASE analysis result
# check variants with more than one max rows
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# be careful of the binom.test result for variants with more than one max rows, because this can cause non-integer readcount 
# in the above calculation which can make errors for binom.test
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tmp<-allele_count1[!is.na(allele_count1$x3_max_num) & !is.na(allele_count1$sum_max_num),]
tmp[tmp$x3_max_num>1|tmp$sum_max_num>1,]
#compare result of two basis: sum of ref readcount and alternate readcount, sequencing depth
table(tmp$Alternate_readcounts_sum_max==tmp$Alternate_readcounts_x3_max)
table(tmp$Reference_readcounts_sum_max==tmp$Reference_readcounts_x3_max);tmp[tmp$Reference_readcounts_sum_max!=tmp$Reference_readcounts_x3_max,]
table(tmp$pval_sum_max==tmp$pval_x3_max);tmp[tmp$pval_sum_max!=tmp$pval_x3_max,]
#p value process: assign NA to p value of rows without enough read count; ajusted with BH
tmpind<-(allele_count1$Reference_readcounts_x3_max+allele_count1$Alternate_readcounts_x3_max)<6; tmpind[is.na(tmpind)]<-F
allele_count1[tmpind,c('pval_sum_max','pval_x3_max')]<-NA; allele_count1$pval_fdr<-p.adjust(allele_count1$pval_x3_max,method = 'BH')
#classification of variant ASE status
classification<-rep('Insufficient read count data',nrow(allele_count1))
classification[(!is.na(allele_count1$pval_fdr))&(allele_count1$pval_fdr>=0.15)]<-'Not ASE'
classification[(!is.na(allele_count1$pval_fdr))&(allele_count1$pval_fdr<0.15)&(allele_count1$pval_fdr>=0.05)]<-'Suggestive'
classification[(!is.na(allele_count1$pval_fdr))&(allele_count1$pval_fdr<0.05)]<-'Significant'
allele_count1$class<-classification
write.table(allele_count1[order(allele_count1$pval_fdr),],file.path(outdir,'ASE_basic-allele_count.txt'),row.names = F,col.names = T,quote = F,sep = '\t')

### plot of ASE analysis result
# barplot for count of sample-variant across different classes
plot_data1<-data.frame(rep('Count_All_Read',4),c('Insufficient read count data','Not ASE','Suggestive','Significant'),c(nrow(variant)-sum(allele_count1$class%in%c('Not ASE','Suggestive','Significant')),sum(allele_count1$class=='Not ASE'),sum(allele_count1$class=='Suggestive'),sum(allele_count1$class=='Significant')),stringsAsFactors = F)
colnames(plot_data1)<-c('Class','Classification','Countss'); plot_data1$Classification<-factor(plot_data1$Classification,levels = plot_data1$Classification)
plot_data1$Ratio<-plot_data1$Countss/sum(plot_data1$Countss)
for (i in 1:nrow(plot_data1)){plot_data1[i,'labloc']<-sum(plot_data1$Countss)-(sum(plot_data1[1:i,'Countss'])-plot_data1[i,'Countss']/2)}
writexl::write_xlsx(plot_data1,file.path(outdir,'ASE_basic-Sample_variant_number_of_ASE.xlsx'),format_headers = F)
ASE_status_color<-c('Significant'=rgb(97,89,164,maxColorValue = 255),'Suggestive'=rgb(206,74,8,maxColorValue = 255),'Not ASE'=rgb(29,143,100,maxColorValue = 255),'Insufficient read count data'='gray70')

p<-ggplot(plot_data1,aes(x=Class,y=Countss,fill=Classification))+geom_bar(stat = 'identity')
p<-p+geom_text(aes(y=labloc,label=paste(c('Insufficient read count data','Not ASE','Suggestive','Significant'),Countss,sep = ',')),size=12)
p<-p+scale_fill_manual(values = ASE_status_color)
p<-p+theme_minimal()+theme(axis.title = element_blank(),axis.text = element_blank(),panel.grid = element_blank(),legend.position = 'none')
p
ggsave(file.path(outdir,'ASE_basic-Sample_variant_number_of_ASE.pdf'),h=7,w=9,useDingbat=F)
#plot of reference read count vs alternate read count for all sample-variant
p<-ggplot(allele_count1,aes(x=Reference_readcounts_x3_max,y=Alternate_readcounts_x3_max))
p<-p+geom_point(aes(color=class),alpha=0.5,stroke=0,size=4)
p<-p+scale_color_manual(name='ASE Status',values = ASE_status_color)
p<-p+labs(x='Reference Read Count',y='Alternate Read Count')
p<-p+theme_bw()+theme(axis.text = element_text(size = 17),axis.title = element_text(size = 20),
                      legend.text = element_text(size = 17),legend.title = element_text(size = 20),
                      plot.margin = unit(c(0,0,0,0),'cm'))
p
ggsave(file.path(outdir,'ASE_basic-Read_Count_of_Different_ASE_Classes_for_All_Variant.pdf'),h=7,w=9,useDingbat=F)
p<-ggplot(allele_count1,aes(x=log10(Reference_readcounts_x3_max),y=log10(Alternate_readcounts_x3_max)))
p<-p+geom_point(aes(color=class),alpha=0.5,stroke=0,size=4)
p<-p+scale_color_manual(name='ASE Status',values = ASE_status_color)
p<-p+labs(x='Log10(Reference Read Count)',y='Log10(Alternate Read Count)')
p<-p+scale_x_continuous(limits = c(-1,max(log10(allele_count1$Reference_readcounts_x3_max))))+scale_y_continuous(limits = c(-1,max(log10(allele_count1$Alternate_readcounts_x3_max))))
p<-p+theme_bw()+theme(axis.text = element_text(size = 17),axis.title = element_text(size = 20),
                      legend.text = element_text(size = 17),legend.title = element_text(size = 20),
                      plot.margin = unit(c(0,0,0,0),'cm'))
p
ggsave(file.path(outdir,'ASE_basic-Log_Read_Count_of_Different_ASE_Classes_for_All_Variant.pdf'),h=7,w=9,useDingbat=F)
#plot of reference read count vs alternate read count for sample-variant with enough readcount
allele_count1_enough_read<-allele_count1[!is.na(allele_count1$pval_fdr),]
p<-ggplot(allele_count1_enough_read,aes(x=Reference_readcounts_x3_max,y=Alternate_readcounts_x3_max))
p<-p+geom_point(aes(color=class),alpha=0.5,stroke=0,size=4)
p<-p+scale_color_manual(name='ASE Status',values = ASE_status_color)
p<-p+labs(x='Reference Read Count',y='Alternate Read Count')
p<-p+theme_bw()+theme(axis.text = element_text(size = 17),axis.title = element_text(size = 20),
           legend.text = element_text(size = 17),legend.title = element_text(size = 20),
           plot.margin = unit(c(0,0,0,0),'cm'))
p
ggsave(file.path(outdir,'ASE_basic-Read_Count_of_Different_ASE_Classes_for_Variant_with_Enough_Read.pdf'),h=7,w=9,useDingbat=F)
p<-ggplot(allele_count1_enough_read,aes(x=log10(Reference_readcounts_x3_max),y=log10(Alternate_readcounts_x3_max)))
p<-p+geom_point(aes(color=class),alpha=0.5,stroke=0,size=4)
p<-p+scale_color_manual(name='ASE Status',values = ASE_status_color)
p<-p+labs(x='Log10(Reference Read Count)',y='Log10(Alternate Read Count)')
p<-p+scale_x_continuous(limits = c(-1,max(log10(allele_count1_enough_read$Reference_readcounts_x3_max))))+scale_y_continuous(limits = c(-1,max(log10(allele_count1_enough_read$Alternate_readcounts_x3_max))))
p<-p+theme_bw()+theme(axis.text = element_text(size = 17),axis.title = element_text(size = 20),
           legend.text = element_text(size = 17),legend.title = element_text(size = 20),
           plot.margin = unit(c(0,0,0,0),'cm'))
p
ggsave(file.path(outdir,'ASE_basic-Log_Read_Count_of_Different_ASE_Classes_for_Variant_with_Enough_Read.pdf'),h=7,w=9,useDingbat=F)
