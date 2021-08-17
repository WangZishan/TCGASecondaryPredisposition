rm(list = ls()); tmp<-lapply(c("ggplot2", "ggrepel", "scales", "statmod"), library, character.only = TRUE)
setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition')
source("analysis/functions/general_function.R")
indir_ASE<-outdir<-"C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition/analysis/allele_expression/SPV"

### load in related info
# load in the result of ASE analysis
# putting sample-variants belonging to multiple genes into different rows, where each row corresponds to one gene
tmpallele_count1<-read.table(file.path(indir_ASE,"ASE_basic-allele_count.txt"),header = T,sep = "\t",stringsAsFactors = F)
tmpind<-grep(';',tmpallele_count1$Symbol_Merge)
if(length(tmpind)!=0){
  tmpallele_count1_1<-tmpallele_count1[setdiff(1:nrow(tmpallele_count1),tmpind),]; tmpallele_count1_2<-tmpallele_count1[tmpind,]
  tmpallele_count1_2_dup<-matrix(nrow = 0,ncol = ncol(tmpallele_count1_2))
  for (i in 1:nrow(tmpallele_count1_2)) {
    tmp0<-matrix(nrow = 0,ncol = ncol(tmpallele_count1_2)); tmp1<-unlist(strsplit(tmpallele_count1_2[i,'Symbol_Merge'],';'))
    for (j in 1:length(tmp1)) {tmp0<-rbind(tmp0,tmpallele_count1_2[i,])}# for
    tmp0$Symbol_Merge<-tmp1; tmpallele_count1_2_dup<-rbind(tmpallele_count1_2_dup,tmp0)
  }# for
  allele_count1<-rbind(tmpallele_count1_1,tmpallele_count1_2_dup)
}else{allele_count1<-tmpallele_count1}
#variant function class
variant_function<-unique(allele_count1$Function)
#functions used in this analysis
variant_function_transform<-function(x){
  tmp<-rep(0,length(variant_function)); names(tmp)<-variant_function
  x_u<-unique(x); for (z in 1:length(x_u)) {tmp[x_u[z]]<-sum(x==x_u[z])}
  tmp<-tmp[tmp!=0]; return(paste(paste(names(tmp),tmp,sep = ","),collapse = "|"))
}#function

### Gene enrichment analysis of ASE variant using fisher test
genes<-unique(allele_count1$Symbol_Merge)
m<-as.data.frame(matrix(nrow = length(genes),ncol = 15, dimnames = list(NULL,c("Symbol_Merge","Variant_num","Variant_Type","Variant_enoughread_num","Variant_enoughread_Type","variant_fdr05_gene_pval","variant_fdr05_gene_pval_power","variant_fdr05_gene_fdr","Variant_fdr05_sig_num","Variant_fdr05_sig_type","variant_fdr15_gene_pval","variant_fdr15_gene_pval_power","variant_fdr15_gene_fdr","Variant_fdr15_sig_num","Variant_fdr15_sig_type"))),stringsAsFactors=F)
for (i in 1:length(genes)) {
  m[i,"Symbol_Merge"]<-genes[i]
  
  tmpallele_count1<-allele_count1[allele_count1$Symbol_Merge==genes[i],]
  m[i,"Variant_num"]<-nrow(tmpallele_count1)
  m[i,"Variant_Type"]<-variant_function_transform(tmpallele_count1$Function)
  m[i,"Variant_enoughread_num"]<-sum(tmpallele_count1$class%in%c("Significant","Suggestive","Not ASE"))
  m[i,"Variant_enoughread_Type"]<-variant_function_transform(tmpallele_count1[tmpallele_count1$class%in%c("Significant","Suggestive","Not ASE"),"Function"])
  m[i,"Variant_fdr05_sig_num"]<-sum(tmpallele_count1$class=="Significant")
  m[i,"Variant_fdr05_sig_type"]<-variant_function_transform(tmpallele_count1[tmpallele_count1$class=="Significant","Function"])
  m[i,"Variant_fdr15_sig_num"]<-sum(tmpallele_count1$class%in%c("Significant","Suggestive"))
  m[i,"Variant_fdr15_sig_type"]<-variant_function_transform(tmpallele_count1[tmpallele_count1$class%in%c("Significant","Suggestive"),"Function"])
  
  tmpm<-matrix(nrow = 2,ncol = 2,dimnames = list(c("In_gene","Out_gene"),c("ASE","ASEnot")))
  tmpm[1,1]<-sum(tmpallele_count1$class=="Significant");  tmpm[1,2]<-sum(tmpallele_count1$class%in%c("Not ASE","Suggestive"))
  tmpm[2,1]<-sum(allele_count1$class=="Significant")-tmpm[1,1];  tmpm[2,2]<-sum(allele_count1$class%in%c("Not ASE","Suggestive"))-tmpm[1,2]
  if(sum(tmpm[1,])!=0){
    m[i,"variant_fdr05_gene_pval"]<-fisher.test(tmpm)$p.value
    m[i,"variant_fdr05_gene_pval_power"]<-power.fisher.test(tmpm[1,1]/sum(tmpm[1,]),tmpm[2,1]/sum(tmpm[2,]),sum(tmpm[1,]),sum(tmpm[2,]),nsim = 1000)
  }#for if
  
  tmpm1<-matrix(nrow = 2,ncol = 2,dimnames = list(c("In_gene","Out_gene"),c("ASE","ASEnot")))
  tmpm1[1,1]<-sum(tmpallele_count1$class%in%c("Significant","Suggestive"));  tmpm1[1,2]<-sum(tmpallele_count1$class=="Not ASE")
  tmpm1[2,1]<-sum(allele_count1$class%in%c("Significant","Suggestive"))-tmpm1[1,1];  tmpm1[2,2]<-sum(allele_count1$class=="Not ASE")-tmpm1[1,2]
  if(sum(tmpm1[1,])!=0){
    m[i,"variant_fdr15_gene_pval"]<-fisher.test(tmpm1)$p.value
    m[i,"variant_fdr15_gene_pval_power"]<-power.fisher.test(tmpm1[1,1]/sum(tmpm1[1,]),tmpm1[2,1]/sum(tmpm1[2,]),sum(tmpm1[1,]),sum(tmpm1[2,]),nsim = 1000)
  }#for if
  
  if(i==1){cat(paste(length(genes),i,sep = " "))}else{cat(paste(" ",i,sep=""))}
}#for i

### process gene enrichment analysis result
#assign NA to pvalue of variants without enough read
m$Variant_fdr05_sig_ratio<-m$Variant_fdr05_sig_num/m$Variant_enoughread_num; m$Variant_fdr15_sig_ratio<-m$Variant_fdr15_sig_num/m$Variant_enoughread_num
tmppval<-m$variant_fdr05_gene_pval; tmppval[!(m$Variant_fdr05_sig_num>=3 & (m$Variant_fdr05_sig_ratio>0.7))]<-NA
m$variant_fdr05_gene_fdr<-p.adjust(tmppval,method = "BH")
tmppval<-m$variant_fdr15_gene_pval; tmppval[!(m$Variant_fdr15_sig_num>=3 & (m$Variant_fdr15_sig_ratio>0.7))]<-NA
m$variant_fdr15_gene_fdr<-p.adjust(tmppval,method = "BH")
m<-m[order(m$variant_fdr05_gene_fdr),]
# classification of gene enrichment
m$Classification<-"None"
m$Classification[!is.na(m$variant_fdr05_gene_fdr) & m$variant_fdr05_gene_fdr<0.05]<-"Significant"
m$Classification[!is.na(m$variant_fdr05_gene_fdr) & m$variant_fdr05_gene_fdr>=0.05 & m$variant_fdr05_gene_fdr<0.15]<-"Suggestive"
write.table(m,file.path(outdir,"ASE_Gene_Enrichment.txt"),row.names = F,col.names = T,quote = F,sep = "\t")
# Count the number of genes enriched with ASE sample-variants
m_sta<-matrix(nrow = 6,ncol = 2,dimnames = list(c("variant_fdr05_gene_pval05","variant_fdr05_gene_fdr05","variant_fdr05_gene_fdr15","variant_fdr15_gene_pval05","variant_fdr15_gene_fdr05","variant_fdr15_gene_fdr15"),c("ASE_Variant_Enriched_Gene_number","Not_ASE_Variant_Enriched_Gene_Number")))
m_sta[1,1]<-sum_notna(m$variant_fdr05_gene_pval<0.05);m_sta[2,1]<-sum_notna(m$variant_fdr05_gene_fdr<0.05);m_sta[3,1]<-sum_notna(m$variant_fdr05_gene_fdr<0.15)
m_sta[4,1]<-sum_notna(m$variant_fdr15_gene_pval<0.05);m_sta[5,1]<-sum_notna(m$variant_fdr15_gene_fdr<0.05);m_sta[6,1]<-sum_notna(m$variant_fdr15_gene_fdr<0.15)
m_sta[,2]<-length(genes)-m_sta[,1]
write.table(m_sta,file.path(outdir,"ASE_Gene_Enrichment-sta.txt"),row.names = T,col.names = T,quote = F,sep = "\t")

### plot of gene enrichment analysis result
# plot for global view of gene enrichment analysis
m$size[m$Classification=="Significant"]<--log2(m$variant_fdr05_gene_fdr[m$Classification=="Significant"])*2
m$size[m$Classification=="Suggestive"]<-4; m$size[m$Classification=="None"]<-0.5
p<-ggplot(m,aes(x=Variant_enoughread_num,y=Variant_fdr05_sig_ratio,color=Classification))+geom_point(aes(size=size))
p<-p+geom_label_repel(aes(label=ifelse(Classification=="None",NA,Symbol_Merge),size=35),show.legend = F)
p<-p+labs(x="Number of varinats with enough readcount\nlocated at the gene region", y="Proportion of ASE variant",color="ASE SPV\nEnrichment",size="Enrichment level\n[-log(FDR)]")
p<-p+theme_bw()
p<-p+theme(axis.text = element_text(size = 13),axis.title = element_text(size = 16),
           legend.text = element_text(size = 13),legend.title = element_text(size=16),
           legend.position = c(0.82,0.7))
p<-p+scale_color_manual(values = c("Significant"="#F8766D","Suggestive"="#00BA38","None"="#619CFF"))
p
ggsave(h=6,w=6,useDingbat=F,filename = file.path(outdir,"ASE_Gene_Enrichment-global_view.pdf"))
# count variant type for genes enriched with ASE
msig_ori<-rbind(m[m$Classification=="Significant",],m[m$Classification=="Suggestive",])
variant_type<-setdiff(unique(unlist(strsplit(unlist(strsplit(msig_ori$Variant_Type,"\\|")),","))),as.character(1:10000))
msig<-data.frame(expand.grid(variant_type,msig_ori$Symbol_Merge,stringsAsFactors = F),0,0,0,stringsAsFactors = F)
colnames(msig)<-c("variant_type","gene","sig","not_sig","all")
for (i in 1:nrow(msig)) {
  tmptmp<-msig_ori[msig_ori$Symbol_Merge==msig[i,"gene"],]
  tmpind<-unlist(gregexpr(msig[i,"variant_type"],tmptmp[,"Variant_Type"]))
  if(tmpind!=-1){
    tmptmpv<-tmptmp[,"Variant_Type"]; tmptmpv1<-unlist(strsplit(tmptmpv,"\\|")); tmptmpv2<-tmptmpv1[grep(msig[i,"variant_type"],tmptmpv1)]
    if(length(tmptmpv2)!=0){msig[i,"all"]<-as.numeric(unlist(strsplit(tmptmpv2,","))[2])}#for if
    tmptmpv<-tmptmp[,"Variant_fdr05_sig_type"]; tmptmpv1<-unlist(strsplit(tmptmpv,"\\|")); tmptmpv2<-tmptmpv1[grep(msig[i,"variant_type"],tmptmpv1)]
    if(length(tmptmpv2)!=0){msig[i,"sig"]<-as.numeric(unlist(strsplit(tmptmpv2,","))[2])}#for if
  }#for if
}#for i
msig[,"not_sig"]<-msig[,"all"]-msig[,"sig"]
# plot for variant type count
tmpplot_data1<-msig[,c("variant_type","gene","sig")]; tmpplot_data2<-msig[,c("variant_type","gene","not_sig")]
colnames(tmpplot_data1)<-colnames(tmpplot_data2)<-c("variant_type","gene","freq")
tmpplot_data2$freq<-tmpplot_data2$freq*(-1); plot_data<-rbind(tmpplot_data1,tmpplot_data2)
p<-ggplot(plot_data,aes(x=gene,y=freq,fill=variant_type))+geom_bar(stat = "identity")+geom_hline(yintercept = 0)
p<-p+theme_minimal()+labs(x="Genes enriched with significant ASE variants",y="Count of not significant or significant variants")
p<-p+scale_y_continuous(breaks = seq(min(plot_data$freq),max(plot_data$freq)+1,2))
p<-p+theme(axis.text.x = element_text(size = 17,angle = 45,vjust = 1,hjust = 1),axis.text.y = element_text(size = 15),axis.title = element_text(size = 20),
           panel.grid.major.y = element_line(color = "grey"),panel.grid.major.x = element_blank(),panel.grid.minor.y = element_blank(),
           plot.margin = unit(c(0,0,0,0),"cm"),legend.text = element_text(size = 15),legend.title = element_text(size = 18))
p
ggsave(h=7,w=8,useDingbat=F,filename = file.path(outdir,"ASE_Gene_Enrichment-different_type_count.pdf"))
# plot for ASE status count of genes enriched with ASE
allele_count_selected = allele_count1[(allele_count1$Symbol_Merge %in% msig_ori$Symbol_Merge) & (allele_count1$class %in% c("Not ASE","Suggestive","Significant")),]
allele_count_selected$Symbol_Merge = factor(allele_count_selected$Symbol_Merge,levels = msig_ori$Symbol_Merge)
allele_count_selected$class<-factor(allele_count_selected$class,levels = c("Not ASE","Suggestive","Significant"))
ASE_status_color<-c("Significant"=rgb(97,89,164,maxColorValue = 255),"Suggestive"=rgb(206,74,8,maxColorValue = 255),"Not ASE"=rgb(29,143,100,maxColorValue = 255))
p<-ggplot(allele_count_selected,aes(x=Symbol_Merge,fill=class))+geom_bar(stat = "count")
p<-p+theme_bw()
p<-p+labs(x="Genes",y="Variant Counts")+scale_y_continuous(breaks = seq(0,max(table(allele_count_selected$Symbol_Merge))+1,2))
p<-p+theme(axis.text.x = element_text(size = 18,angle = 90,hjust = 1,vjust = 0.5),axis.text.y = element_text(size = 18),
           axis.title = element_text(size = 20),
           legend.title = element_text(size = 18),legend.text = element_text(size = 15),
           panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank())
p<-p+scale_fill_manual(name="ASE Status",values = ASE_status_color)
p
ggsave(h=6,w=9,useDingbat=F,filename = file.path(outdir,"ASE_Gene_Enrichment-ase_status_count.pdf"))

