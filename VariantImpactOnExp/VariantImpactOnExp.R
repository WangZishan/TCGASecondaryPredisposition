rm(list = ls()); setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition')
source('analysis/functions/stat_functions.R')
outdir<-'C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition/analysis/VariantImpactOnExp/exp_percentile'

### Load in related information
# tumor.purity
tp <- as.data.frame(readxl::read_xlsx('analysis/data/TumorPurity/41467_2015_BFncomms9971_MOESM1236_ESM_d_d.xlsx'),stringsAsFactors=F)
vars <- c('ESTIMATE','ABSOLUTE','LUMP','IHC','CPE'); for (i in 1:length(vars)) {tp[,vars[i]] <- as.numeric(tp[,vars[i]])}# for
# primary component information
pc<-read.table('analysis/data/primary_component_data_preprocess/primary_component.txt',header = T,row.names = 1,sep = '\t')[,c('PC1','PC2')]
# clinical information
clin_complete = read.table('../Huang_lab/Huang_lab_data/TCGA_PanCanAtlas_2018/GDC_Data_2018/clinical_PANCAN_patient_with_followup.tsv', header=T, quote = '', sep='\t', fill =T, stringsAsFactors=FALSE)
clin_complete<-clin_complete[!clin_complete$age_at_initial_pathologic_diagnosis=='[Not Available]',]
clin_complete$age_at_initial_pathologic_diagnosis<-as.numeric(clin_complete$age_at_initial_pathologic_diagnosis)
clin_brief = clin_complete[,c('bcr_patient_barcode','acronym','age_at_initial_pathologic_diagnosis','gender')]
# read in expression profile
exp_normal3<-read.table('analysis/VariantImpactOnExp/data/expression_of_variant_genes_preprocess/normal_expression_of_variant_genes-pc-clinical-ecdf.txt',row.names = 1,header = T,sep = '\t',stringsAsFactors = F)
exp_cancer3<-read.table('analysis/VariantImpactOnExp/data/expression_of_variant_genes_preprocess/cancer_expression_of_variant_genes-pc-clinical-ecdf.txt',row.names = 1,header = T,sep = '\t',stringsAsFactors = F)
colnames(exp_normal3)<-gsub('\\.','-',colnames(exp_normal3)); colnames(exp_cancer3)<-gsub('\\.','-',colnames(exp_cancer3))
exp_normal_sample1<-read.table('analysis/VariantImpactOnExp/data/expression_of_variant_genes_preprocess/normal_expression_of_variant_genes-pc-clinical-sampleinfo.txt',header = T,sep = '\t',stringsAsFactors = F)
exp_cancer_sample1<-read.table('analysis/VariantImpactOnExp/data/expression_of_variant_genes_preprocess/cancer_expression_of_variant_genes-pc-clinical-sampleinfo.txt',header = T,sep = '\t',stringsAsFactors = F)
# variant information
variant<-read.table('analysis/data/variant_data/variants_pca_frq_nonCancer_pathogenic.txt',header = T,sep = '\t',stringsAsFactors = F,quote = '')
variant1<-variant[variant$bcr_patient_barcode%in%colnames(exp_cancer3),]
variant2<-variant[variant$bcr_patient_barcode%in%colnames(exp_normal3),]

### impact of variant on gene expression using linear regression model
#individual cancer in cancer samples
r1<-1
results_list1<-list()
variant_sym_exp1<-unique(variant1$Symbol_Merge)
for (i in 1:length(variant_sym_exp1)) {
  tmpvariant_sym_exp1<-variant_sym_exp1[i]
  tmpvariant1<-variant1[variant1$Symbol_Merge==tmpvariant_sym_exp1,]
  tmpcs<-unique(tmpvariant1$Cancer_Type)
  for (j in 1:length(tmpcs)) {
    tmptmpcs<-tmpcs[j]
    tmptmpvariant1<-tmpvariant1[tmpvariant1$Cancer_Type==tmptmpcs,]
    
    tmpexp<-t(exp_cancer3[tmpvariant_sym_exp1,exp_cancer_sample1[,2]==tmptmpcs])[,1]
    if(sum(is.na(tmpexp))==0){
      tmpsample<-names(tmpexp)

      tmpcarrier<-rep(F,length(tmpexp))
      tmpcarrier[tmpsample%in%tmptmpvariant1$bcr_patient_barcode]<-T
      
      tmpgender_age<-clin_brief[match(tmpsample,clin_brief$bcr_patient_barcode),c('gender','age_at_initial_pathologic_diagnosis')]
      tmppc<-pc[tmpsample,]
      
      tmptp <- tp[match(tmpsample,tp$bcr_patient_barcode),'CPE']
      
      tmpdata<-data.frame(tmpexp,tmpcarrier,tmpgender_age,tmppc,tmptp,stringsAsFactors = F)
      colnames(tmpdata)<-c('exp','mutation','gender','age','PC1','PC2','TumorPurity')
      for (k in c(1,4,5,6,7)) {tmpdata[,k]<-as.numeric(tmpdata[,k])}#for k
      tmpdata[,3]<-as.character(tmpdata[,3])
      
      fit_linear<-run_glm(data=tmpdata,yi='exp',xi='mutation',ytype = 'Continuous',covi = c('age','gender','PC1','PC2','TumorPurity'))
      cancer_stat = data.frame(cbind(tmpvariant_sym_exp1, tmptmpcs, length(tmpcarrier), sum(tmpcarrier), fit_linear))
      results_list1[[r1]] = cancer_stat
      r1=r1+1
    }#for if
  }#for j
}#for i

#individual cancer in normal samples
tmp<-as.data.frame(table(exp_normal_sample1[,2]),stringsAsFactors = F)
cs_not<-tmp[tmp[,2]<6,1]
r2<-1
results_list2<-list()
variant_sym_exp2<-unique(variant2$Symbol_Merge)
for (i in 1:length(variant_sym_exp2)) {
  tmpvariant_sym_exp2<-variant_sym_exp2[i]
  tmpvariant2<-variant2[variant2$Symbol_Merge==tmpvariant_sym_exp2,]
  tmpcs<-unique(tmpvariant2$Cancer_Type)
  tmpcs<-setdiff(tmpcs,cs_not)
  tmpcs_l<-length(tmpcs)
  if(tmpcs_l!=0){
    for (j in 1:tmpcs_l) {
      tmptmpcs<-tmpcs[j]
      tmptmpvariant2<-tmpvariant2[tmpvariant2$Cancer_Type==tmptmpcs,]
      
      tmpexp<-t(exp_normal3[tmpvariant_sym_exp2,exp_normal_sample1[,2]==tmptmpcs])[,1]
      if(sum(is.na(tmpexp))==0){
        tmpsample<-names(tmpexp)

        tmpcarrier<-rep(F,length(tmpexp))
        tmpcarrier[tmpsample%in%tmptmpvariant2$bcr_patient_barcode]<-T
        
        tmpgender_age<-clin_brief[match(tmpsample,clin_brief$bcr_patient_barcode),c('gender','age_at_initial_pathologic_diagnosis')]
        tmppc<-pc[tmpsample,]
        
        tmptp <- tp[match(tmpsample,tp$bcr_patient_barcode),'CPE']
        
        tmpdata<-data.frame(tmpexp,tmpcarrier,tmpgender_age,tmppc,tmptp,stringsAsFactors = F)
        colnames(tmpdata)<-c('exp','mutation','gender','age','PC1','PC2','TumorPurity')
        for (k in c(1,4,5,6,7)) {tmpdata[,k]<-as.numeric(tmpdata[,k])}#for k
        tmpdata[,3]<-as.character(tmpdata[,3])
        
        fit_linear<-run_glm(data=tmpdata,yi='exp',xi='mutation',ytype = 'Continuous',covi = c('age','gender','PC1','PC2','TumorPurity'))
        cancer_stat = data.frame(cbind(tmpvariant_sym_exp2, tmptmpcs, length(tmpcarrier), sum(tmpcarrier), fit_linear))
        results_list2[[r2]] = cancer_stat
        r2=r2+1
      }#for if
    }#for j
  }#for if
}#for i

#pan cancer in cancer samples
r3<-1
results_list3<-list()
variant_sym_exp1<-unique(variant1$Symbol_Merge)

tmpsample<-exp_cancer_sample1$patient_name
tmpgender_age<-clin_brief[match(tmpsample,clin_brief$bcr_patient_barcode),c('gender','age_at_initial_pathologic_diagnosis')]
tmppc<-pc[tmpsample,]
tmpcancertyp<-exp_cancer_sample1$cancer_name
tmptp <- tp[match(tmpsample,tp$bcr_patient_barcode),'CPE']
tmptmpdata<-data.frame(tmpgender_age,tmppc,tmpcancertyp,tmptp,stringsAsFactors = F)
colnames(tmptmpdata)<-c('gender','age','PC1','PC2','cancertype','TumorPurity')
for (i in 1:length(variant_sym_exp1)) {
  tmpvariant_sym_exp1<-variant_sym_exp1[i]
  tmpvariant1<-variant1[variant1$Symbol_Merge==tmpvariant_sym_exp1,]
  
  if(tmpvariant_sym_exp1%in%rownames(exp_cancer3)){
    tmpexp<-t(exp_cancer3[tmpvariant_sym_exp1,])[,1]
    tmpsample<-names(tmpexp)
    
    tmpcarrier<-rep(F,length(tmpexp))
    tmpcarrier[tmpsample%in%tmpvariant1$bcr_patient_barcode]<-T
    
    tmpdata<-data.frame(tmpexp,tmpcarrier,tmptmpdata,stringsAsFactors = F)
    colnames(tmpdata)<-c('exp','mutation','gender','age','PC1','PC2','cancertype','TumorPurity')
    
    for (k in c(1,4,5,6,8)) {tmpdata[,k]<-as.numeric(tmpdata[,k])}#for k
    for (k in c(3,7)) {tmpdata[,k]<-as.character(tmpdata[,k])}#for k
    
    fit_linear<-run_glm(data=tmpdata,yi='exp',xi='mutation',ytype = 'Continuous',covi = c('age','gender','PC1','PC2','cancertype','TumorPurity'))
    cancer_stat = data.frame(cbind(tmpvariant_sym_exp1, 'pan-cancer', length(tmpcarrier), sum(tmpcarrier), fit_linear))
    results_list3[[r3]] = cancer_stat
    r3=r3+1
  }#for if
}#for i

#pan cancer in normal samples
r4<-1
results_list4<-list()
variant_sym_exp2<-unique(variant2$Symbol_Merge)

tmpsample<-exp_normal_sample1$patient_name
tmpgender_age<-clin_brief[match(tmpsample,clin_brief$bcr_patient_barcode),c('gender','age_at_initial_pathologic_diagnosis')]
tmppc<-pc[tmpsample,]
tmpcancertyp<-exp_normal_sample1$cancer_name
tmptp <- tp[match(tmpsample,tp$bcr_patient_barcode),'CPE']
tmptmpdata<-data.frame(tmpgender_age,tmppc,tmpcancertyp,tmptp,stringsAsFactors = F)
colnames(tmptmpdata)<-c('gender','age','PC1','PC2','cancertype','TumorPurity')
for (i in 1:length(variant_sym_exp2)) {
  tmpvariant_sym_exp2<-variant_sym_exp2[i]
  tmpvariant2<-variant2[variant2$Symbol_Merge==tmpvariant_sym_exp2,]
  
  if(tmpvariant_sym_exp2%in%rownames(exp_normal3)){
    tmpexp<-t(exp_normal3[tmpvariant_sym_exp2,])[,1]
    tmpsample<-names(tmpexp)
    
    tmpcarrier<-rep(F,length(tmpexp))
    tmpcarrier[tmpsample%in%tmpvariant2$bcr_patient_barcode]<-T
    
    tmpdata<-data.frame(tmpexp,tmpcarrier,tmptmpdata,stringsAsFactors = F)
    colnames(tmpdata)<-c('exp','mutation','gender','age','PC1','PC2','cancertype','TumorPurity')
    
    for (k in c(1,4,5,6,8)) {tmpdata[,k]<-as.numeric(tmpdata[,k])}#for k
    for (k in c(3,7)) {tmpdata[,k]<-as.character(tmpdata[,k])}#for k
    
    fit_linear<-run_glm(data=tmpdata,yi='exp',xi='mutation',ytype = 'Continuous',covi = c('age','gender','PC1','PC2','cancertype','TumorPurity'))
    cancer_stat = data.frame(cbind(tmpvariant_sym_exp2, 'pan-cancer', length(tmpcarrier), sum(tmpcarrier), fit_linear))
    results_list4[[r4]] = cancer_stat
    r4=r4+1
  }#for if
}#for i

### Process result and output
tt1 = do.call(rbind,results_list1)
colnames(tt1) = c('gene','cancer','num_patient','num_mutation_patient','y',
                  'y_type','x','degrees_freedom','deviance','residual_degrees_freedom',
                  'residual_deviance','F','P_Chi','coefficient','covariates');
tt1$FDR<-p.adjust(tt1$P_Chi,method = 'BH')
tmp<-tt1$P_Chi
tmp[tt1$num_mutation_patient<3]<-NA
tmp1<-p.adjust(tmp,method = 'BH')
tt1$P_Chi1<-tmp
tt1$FDR1<-tmp1
write.table(tt1,file.path(outdir,'VariantImpactOnExp-individual_cancer.txt'),row.names = F,col.names = T,quote = F,sep = '\t')

tt2 = do.call(rbind,results_list2)
colnames(tt2) = c('gene','cancer','num_patient','num_mutation_patient','y',
                  'y_type','x','degrees_freedom','deviance','residual_degrees_freedom',
                  'residual_deviance','F','P_Chi','coefficient','covariates');
tt2$FDR<-p.adjust(tt2$P_Chi,method = 'BH')
tmp<-tt2$P_Chi
tmp[tt2$num_mutation_patient<3]<-NA
tmp1<-p.adjust(tmp,method = 'BH')
tt2$P_Chi1<-tmp
tt2$FDR1<-tmp1
write.table(tt2,file.path(outdir,'VariantImpactOnExp-individual_normal.txt'),row.names = F,col.names = T,quote = F,sep = '\t')

tt3 = do.call(rbind,results_list3)
colnames(tt3) = c('gene','cancer','num_patient','num_mutation_patient','y',
                  'y_type','x','degrees_freedom','deviance','residual_degrees_freedom',
                  'residual_deviance','F','P_Chi','coefficient','covariates');
tt3$FDR<-p.adjust(tt3$P_Chi,method = 'BH')
tmp<-tt3$P_Chi
tmp[tt3$num_mutation_patient<3]<-NA
tmp1<-p.adjust(tmp,method = 'BH')
tt3$P_Chi1<-tmp
tt3$FDR1<-tmp1
write.table(tt3,file.path(outdir,'VariantImpactOnExp-pan_cancer.txt'),row.names = F,col.names = T,quote = F,sep = '\t')

tt4 = do.call(rbind,results_list4)
colnames(tt4) = c('gene','cancer','num_patient','num_mutation_patient','y',
                  'y_type','x','degrees_freedom','deviance','residual_degrees_freedom',
                  'residual_deviance','F','P_Chi','coefficient','covariates');
tt4$FDR<-p.adjust(tt4$P_Chi,method = 'BH')
tmp<-tt4$P_Chi
tmp[tt4$num_mutation_patient<3]<-NA
tmp1<-p.adjust(tmp,method = 'BH')
tt4$P_Chi1<-tmp
tt4$FDR1<-tmp1
write.table(tt4,file.path(outdir,'VariantImpactOnExp-pan_normal.txt'),row.names = F,col.names = T,quote = F,sep = '\t')

