rm(list = ls()); tmp<-lapply(c('ggplot2','reshape2','grid'), library,character.only=T)

setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition')
outdir<-file.path(getwd(),'analysis','variant_distribution_of_ancestry','variants_pca_frq_nonCancer_pathogenic')
###read in acmg gene list
acmg_gene<-read.table('analysis/data/ACMG_classification2/ACMG_genelist.txt',header = T,sep = '\t',stringsAsFactors = F)[,'ACMG_genelist']
###ancestry order
ancestry<-read.table('analysis/data/ancestry_info/ancestry_patient.txt',header = T,sep = '\t',stringsAsFactors = F)
ancestry_order<-c('Amr','Afr','Eur','Eas','Sas','Mix','Oth')
###read in sample-variant dataset
###putting sample-variant belonging to multiple genes into different rows, where each row corresponds to one gene
tmpvariant<-read.table('analysis/data/variant_data/variants_pca_frq_nonCancer_pathogenic.txt',header = T,sep = '\t',quote = '',stringsAsFactors = F)
tmpind<-grep(';',tmpvariant$Symbol_Merge)
if(length(tmpind)!=0){
  tmpvariant1<-tmpvariant[setdiff(1:nrow(tmpvariant),tmpind),]; tmpvariant2<-tmpvariant[tmpind,]
  tmpvariant2_dup<-matrix(nrow = 0,ncol = ncol(tmpvariant))
  for (i in 1:nrow(tmpvariant2)) {
    tmp0<-matrix(nrow = 0,ncol = ncol(tmpvariant)); tmp1<-unlist(strsplit(tmpvariant2[i,'Symbol_Merge'],';'))
    for (j in 1:length(tmp1)) {tmp0<-rbind(tmp0,tmpvariant2[i,])}# for
    tmp0$Symbol_Merge<-tmp1; tmpvariant2_dup<-rbind(tmpvariant2_dup,tmp0)
  }# for
  variant<-rbind(tmpvariant1,tmpvariant2_dup)
}else{
  variant <- tmpvariant
}

###Count/Percentage of patients with variant in different ancestries
tmp<-unique(variant[,c('bcr_patient_barcode','Symbol_Merge','ancestry')])
count_long<-as.data.frame(table(tmp[,c('Symbol_Merge','ancestry')]),stringsAsFactors = F)
count_wide<-dcast(count_long,formula = Symbol_Merge~ancestry,value.var = 'Freq')
tmp<-t(count_wide[,ancestry_order]); tmp<-t(tmp/as.numeric(table(ancestry$ancestry)[ancestry_order]))
percentage_wide<-data.frame(count_wide$Symbol_Merge,tmp,stringsAsFactors = F); colnames(percentage_wide)<-c('Symbol_Merge',ancestry_order)

percentage_wide_rowsum <- apply(percentage_wide[,setdiff(colnames(percentage_wide),c('Symbol_Merge','Mix','Oth'))], 1, sum)
tmpind <- order(percentage_wide_rowsum,decreasing = T); percentage_wide_rowsum <- percentage_wide_rowsum[tmpind]; percentage_wide <- percentage_wide[tmpind,]
count_wide<-count_wide[match(percentage_wide$Symbol_Merge,count_wide$Symbol_Merge),]
write.table(count_wide,file.path(outdir,'Distribution_for_genes-Count_all.txt'),row.names = F,col.names = T,quote = F,sep = '\t')
write.table(percentage_wide,file.path(outdir,'Distribution_for_genes-Frequency_all.txt'),row.names = F,col.names = T,quote = F,sep = '\t')

###Percentage/Count of patients affected by our variant for a gene subset: top10percent and acmg
#Percentage
percentage_wide_subset1 <- percentage_wide[percentage_wide_rowsum>=percentage_wide_rowsum[length(percentage_wide_rowsum)*0.1],]
percentage_wide_subset1 <- percentage_wide_subset1[percentage_wide_subset1$Symbol_Merge%in%setdiff(percentage_wide_subset1$Symbol_Merge,acmg_gene),]
percentage_wide_subset2 <- percentage_wide[percentage_wide$Symbol_Merge%in%acmg_gene,]
percentage_wide_subset <- rbind(percentage_wide_subset2,percentage_wide_subset1)
write.table(percentage_wide_subset,file.path(outdir,'Distribution_for_genes2-Frequency.txt'),row.names = F,col.names = T,quote = F,sep = '\t')

percentage_long_subset <- melt(percentage_wide_subset,id.vars = 'Symbol_Merge'); colnames(percentage_long_subset) <- c('Symbol_Merge','Ancestry','Frequency')
percentage_long_subset$Symbol_Merge <- factor(percentage_long_subset$Symbol_Merge,levels = rev(percentage_wide_subset$Symbol_Merge))
percentage_long_subset$Ancestry<-factor(percentage_long_subset$Ancestry,levels = ancestry_order)
percentage_long_subset$Frequency <- round(percentage_long_subset$Frequency*100,2)
p<-ggplot(data = percentage_long_subset,aes(x=Ancestry,y=Symbol_Merge,fill=Frequency))+geom_tile(color='grey')
p<-p+theme_minimal()+theme(axis.text.x = element_text(size = 17),axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),axis.title.y = element_blank(),
        panel.grid.major = element_blank(),plot.margin = unit(c(0,0,0,0),'cm'),
        legend.text = element_text(size = 15,angle = 45),legend.title = element_text(size = 18),legend.position = 'bottom')
p<-p+scale_fill_gradient2(low = 'white',high = 'red',name='Frequency [%]')
p<-p+coord_fixed(ratio = 0.25)
p<-p+geom_text(aes(label=ifelse(Frequency!=0,Frequency,NA)),color='black',size=4)
p
ggsave(p,filename = file.path(outdir,'Distribution_for_genes2-Frequency.pdf'),h=7,w=8,useDingbat=F)

#Count
count_wide_subset<-count_wide[match(percentage_wide_subset$Symbol_Merge,count_wide$Symbol_Merge),]
write.table(count_wide_subset,file.path(outdir,'Distribution_for_genes2-Count.txt'),row.names = F,col.names = T,quote = F,sep = '\t')

count_long_subset<-melt(count_wide_subset,id.vars = 'Symbol_Merge'); colnames(count_long_subset)<-c('Symbol_Merge','Ancestry','Count')
count_long_subset$Symbol_Merge<-factor(count_long_subset$Symbol_Merge,levels = rev(count_wide_subset$Symbol_Merge))
count_long_subset$Ancestry<-factor(count_long_subset$Ancestry,levels = ancestry_order)
p<-ggplot(data = count_long_subset,aes(x=Ancestry,y=Symbol_Merge,fill=Count))+geom_tile(color='grey')
p<-p+theme_minimal()+theme(axis.text.x = element_text(vjust = 1,size = 17),axis.text.y = element_text(size = 12),
                           axis.title.x = element_blank(),axis.title.y = element_blank(),
                           panel.grid.major = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
                           plot.margin = unit(c(0,0,0,0),'cm'),legend.text = element_text(size = 15),legend.title = element_text(size = 18),legend.position = 'bottom')
p<-p+scale_fill_gradient2(low = 'white',high = 'red',name='Count')
p<-p+coord_fixed(ratio = 0.25)
p<-p+geom_text(aes(label=ifelse(Count!=0,Count,NA)),color='black',size=4)
p
ggsave(p,filename = file.path(outdir,'Distribution_for_genes2-Count.pdf'),h=7,w=8,useDingbat=F)
