rm(list = ls()); setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition')
tmp<-lapply(c('reshape2','ggplot2','scales'), library, character.only=T)
# directories
indir_gnomad<-'analysis/data/gnomad/bcftools_process/process/'
indir_genes<-'analysis/variant_distribution_of_ancestry/variants_pca_frq_nonCancer_pathogenic/'
outdir<-'analysis/data/gnomad/bcftools_process/plot_ancestry/'

### Load in related information
# Gene list
genes_acmg<-read.table(file.path(indir_genes,'Distribution_for_genes-Count_acmg.txt'),header = T,sep = '\t',stringsAsFactors = F)$Symbol_Merge
genes_top10p<-read.table(file.path(indir_genes,'Distribution_for_genes-Count_top10percent.txt'),header = T,sep = '\t',stringsAsFactors = F)$Symbol_Merge

### plot
# barplot
tmp<-paste(c('non_cancer_AC','non_cancer_AF','non_cancer_AF1'),'_acmg_status',sep = '')
for (i in 1:length(tmp)) {
  tmps<-t(read.table(file.path(indir_gnomad,paste(tmp[i],'.txt',sep = '')),header = T,sep = '\t',stringsAsFactors = F))
  tmps<-tmps[setdiff(1:nrow(tmps),c(grep('Asian',rownames(tmps)),grep('Eur',rownames(tmps)))),]
  tmps_rnames<-rownames(tmps);  tmps1<-data.frame(tmps_rnames,tmps)
  write.table(tmps1,file.path(outdir,paste(tmp[i],'.txt',sep = '')),row.names = F,col.names = T,quote = F,sep = '\t')
  
  tmptmp<-melt(tmps1,id='tmps_rnames')
  tmptmp$tmps_rnames<-factor(tmptmp$tmps_rnames,levels = rev(rownames(tmps)))
  if(length(grep('AC',tmp[i]))!=0){
    for (j in 1:nrow(tmptmp)) {if(tmptmp[j,'variable']=='pathogenic'){tmptmp[j,'loc']<-tmptmp[j,'value']}else{tmptmp[j,'loc']<-sum(tmptmp[tmptmp$tmps_rnames==tmptmp[j,'tmps_rnames'],'value'])}}#for j
    tmptmp1<-tmptmp[tmptmp$variable!='pathogenic',]; tmptmp2<-tmptmp[tmptmp$variable=='pathogenic',]
    p<-ggplot()+geom_bar(data=tmptmp1,aes(x=tmps_rnames,y=loc),fill='#F8766D',stat = 'identity')
    p<-p+geom_bar(data = tmptmp2,aes(x=tmps_rnames,y=loc),fill='#00BFC4',stat = 'identity')
    p<-p+geom_text(data = tmptmp1,aes(x=tmps_rnames,y=loc*0.7,label=value,size=10))
    p<-p+geom_text(data = tmptmp2,aes(x=tmps_rnames,y=loc*0.1,label=value,size=10))
    p<-p+theme_bw()
    p<-p+theme(axis.text.x = element_text(vjust = 1,hjust = 0.5,size = 14),
               axis.text.y = element_text(vjust = 0.5,hjust = 1,size = 14),axis.title = element_blank(),
               legend.position = 'none',panel.grid.major.y = element_blank())
    p<-p+coord_flip()
    p<-p+scale_y_continuous(trans = 'log10')
    p
  }else{
    tmptmp$value<-round(tmptmp$value*100,2)
    for (j in 1:nrow(tmptmp)) {if(tmptmp[j,'variable']=='pathogenic'){tmptmp[j,'loc']<-tmptmp[j,'value']}else{tmptmp[j,'loc']<-sum(tmptmp[tmptmp$tmps_rnames==tmptmp[j,'tmps_rnames'],'value'])}}#for j
    tmptmp1<-tmptmp[tmptmp$variable!='pathogenic',]; tmptmp2<-tmptmp[tmptmp$variable=='pathogenic',]
    p<-ggplot()+geom_bar(data=tmptmp1,aes(x=tmps_rnames,y=loc),fill='#F8766D',stat = 'identity')
    p<-p+geom_bar(data=tmptmp2,aes(x=tmps_rnames,y=loc),fill='#00BFC4',stat = 'identity')
    p<-p+geom_text(data = tmptmp1, aes(x=tmps_rnames,y=loc*0.9,label=value,size=10))
    p<-p+geom_text(data = tmptmp2,aes(x=tmps_rnames,y=loc*0.1,label=value,size=10))
    p<-p+theme_bw()
    p<-p+theme(axis.text.x = element_text(vjust = 1,hjust = 0.5,size=14),
               axis.text.y = element_text(vjust = 0.5,hjust = 1,size = 14),axis.title = element_blank(),
               legend.position = 'none',panel.grid.major.y = element_blank())
    p<-p+coord_flip()
    if(length(grep('AF',tmp[i]))!=0){p<-p+scale_y_reverse()}
    p
  }
  ggsave(file.path(outdir,paste(tmp[i],'.pdf',sep = '')),h=5,w=5)
}#for i

# heatmap plot
options(scipen=999)
tmp<-paste(c('non_cancer_AC','non_cancer_AF','non_cancer_AF1'),'gene',sep = '_')
for (i in 1:length(tmp)) {
  tmps<-read.table(file.path(indir_gnomad,paste(tmp[i],'.txt',sep = '')),header = T,sep = '\t',stringsAsFactors = F)
  tmps<-tmps[,setdiff(1:ncol(tmps),c(grep('Asian',colnames(tmps)),grep('Eur',colnames(tmps))))]
  tmps_rnames<-rownames(tmps);  tmps1<-data.frame(tmps_rnames,tmps,stringsAsFactors = F)
  
  tmps2<-tmps1[as.character(tmps1$tmps_rnames)%in%genes_acmg,]
  tmps3_genes<-setdiff(genes_acmg,rownames(tmps2))
  tmps3<-data.frame(tmps3_genes,matrix(as.integer(0),ncol = ncol(tmps2)-1,nrow = length(tmps3_genes)),stringsAsFactors = F)
  colnames(tmps3)<-colnames(tmps2);  rownames(tmps3)<-tmps3_genes
  tmps4<-rbind(tmps3,tmps2)
  tmptmp<-melt(tmps4,id='tmps_rnames')
  tmptmp$tmps_rnames<-factor(tmptmp$tmps_rnames,levels = rev(genes_acmg))
  if(length(grep('AF',tmp[i]))!=0){tmptmp$value<-round(tmptmp$value,6)*100}#for if 
  p<-ggplot(tmptmp,aes(y=tmps_rnames,x=variable,fill=value))+geom_tile(color='grey',height=1)
  p<-p+geom_text(aes(label=ifelse(value!=0,value,NA),size=10))
  p<-p+theme_minimal()
  p<-p+theme(axis.text.x = element_text(size = 17,angle = 90,hjust = 1,vjust = 0.5),
             axis.text.y = element_text(size = 17),axis.title = element_blank(),legend.position = 'top',
             legend.text = element_text(size = 17,angle = 45,hjust = 1,vjust = 1),legend.title = element_blank(),
             panel.grid = element_blank())
  p<-p+scale_fill_gradient2(low = 'white',high = 'red')
  p
  ggsave(file.path(outdir,paste(tmp[i],'_acmg','.pdf',sep = '')),h=6,w=8)
  write.table(tmps1[genes_acmg,],file.path(outdir,paste(tmp[i],'_acmg','.txt',sep = '')),row.names = F,col.names = T,quote = F,sep = '\t')
  
  tmps2<-tmps1[as.character(tmps1$tmps_rnames)%in%genes_top10p,]
  tmps3_genes<-setdiff(genes_top10p,rownames(tmps2))
  tmps3<-data.frame(tmps3_genes,matrix(as.integer(0),ncol = ncol(tmps2)-1,nrow = length(tmps3_genes)),stringsAsFactors = F)
  colnames(tmps3)<-colnames(tmps2);  rownames(tmps3)<-tmps3_genes
  tmps4<-rbind(tmps3,tmps2)
  tmptmp<-melt(tmps4,id='tmps_rnames')
  tmptmp$tmps_rnames<-factor(tmptmp$tmps_rnames,levels = rev(genes_top10p))
  if(length(grep('AF',tmp[i]))!=0){tmptmp$value<-round(tmptmp$value,6)*100}#for if 
  p<-ggplot(tmptmp,aes(y=tmps_rnames,x=variable,fill=value))+geom_tile(color='grey',height=1)
  p<-p+geom_text(aes(label=ifelse(value!=0,value,NA),size=10))
  p<-p+theme_minimal()
  p<-p+theme(axis.text.x = element_text(size = 17,angle = 90,hjust = 1,vjust = 0.5), 
             axis.text.y = element_text(size = 17),axis.title = element_blank(),legend.position = 'top',
             legend.text = element_text(size = 17,angle = 45,hjust = 1,vjust = 1),legend.title = element_blank(),
             panel.grid = element_blank())
  p<-p+scale_fill_gradient2(low = 'white',high = 'red')
  p
  ggsave(file.path(outdir,paste(tmp[i],'_top10percent','.pdf',sep = '')),h=16,w=8)
  write.table(tmps4[genes_top10p,],file.path(outdir,paste(tmp[i],'_top10percent','.txt',sep = '')),row.names = F,col.names = T,quote = F,sep = '\t')
}#for i

# heatmap plot
options(scipen=999)
tmp<-paste(c('non_cancer_AC','non_cancer_AF','non_cancer_AF1'),'gene',sep = '_')
for (i in 1:length(tmp)) {
  tmps<-read.table(file.path(indir_gnomad,paste(tmp[i],'.txt',sep = '')),header = T,sep = '\t',stringsAsFactors = F)
  tmps<-tmps[,setdiff(1:ncol(tmps),c(grep('Asian',colnames(tmps)),grep('Eur',colnames(tmps))))]
  tmps_rnames<-rownames(tmps);  tmps1<-data.frame(tmps_rnames,tmps,stringsAsFactors = F)
  
  tmps2<-tmps1[as.character(tmps1$tmps_rnames)%in%genes_top10p,]
  tmps3_genes<-setdiff(genes_top10p,rownames(tmps2))
  tmps3<-data.frame(tmps3_genes,matrix(as.integer(0),ncol = ncol(tmps2)-1,nrow = length(tmps3_genes)),stringsAsFactors = F)
  colnames(tmps3)<-colnames(tmps2);  rownames(tmps3)<-tmps3_genes
  tmps4<-rbind(tmps3,tmps2)
  tmpplot1 <- tmps4
  
  tmps2<-tmps1[as.character(tmps1$tmps_rnames)%in%genes_acmg,]
  tmps3_genes<-setdiff(genes_acmg,rownames(tmps2))
  tmps3<-data.frame(tmps3_genes,matrix(as.integer(0),ncol = ncol(tmps2)-1,nrow = length(tmps3_genes)),stringsAsFactors = F)
  colnames(tmps3)<-colnames(tmps2);  rownames(tmps3)<-tmps3_genes
  tmps4<-rbind(tmps3,tmps2)
  tmps4 <- tmps4[setdiff(tmps4$tmps_rnames,tmpplot1$tmps_rnames),]
  tmpplot2 <- tmps4
  tmpplot <- rbind(tmpplot1,tmpplot2)
  
  tmptmp<-melt(tmpplot,id='tmps_rnames')
  tmptmp$tmps_rnames<-factor(tmptmp$tmps_rnames,levels = rev(c(genes_top10p,setdiff(genes_acmg,genes_top10p))))
  if(length(grep('AF',tmp[i]))!=0){tmptmp$value<-round(tmptmp$value,6)*100}#for if 
  p<-ggplot(tmptmp,aes(y=tmps_rnames,x=variable,fill=value))+geom_tile(color='grey',height=1)
  p<-p+geom_text(aes(label=ifelse(value!=0,value,NA),size=10))
  p<-p+theme_minimal()
  p<-p+theme(axis.text.x = element_text(size = 17,angle = 90,hjust = 1,vjust = 0.5),
             axis.text.y = element_text(size = 17),axis.title = element_blank(),legend.position = 'top',
             legend.text = element_text(size = 17,angle = 45,hjust = 1,vjust = 1),legend.title = element_blank(),
             panel.grid = element_blank())
  p<-p+scale_fill_gradient2(low = 'white',high = 'red')
  p
  ggsave(file.path(outdir,paste(tmp[i],'1.pdf',sep = '')),h=14,w=8)
  write.table(tmps1[genes_acmg,],file.path(outdir,paste(tmp[i],'1.txt',sep = '')),row.names = F,col.names = T,quote = F,sep = '\t')
}#for i

# heatmap plot
options(scipen=999)
tmp<-paste(c('non_cancer_AC','non_cancer_AF','non_cancer_AF1'),'gene',sep = '_')
for (i in 1:length(tmp)) {
  tmps<-read.table(file.path(indir_gnomad,paste(tmp[i],'.txt',sep = '')),header = T,sep = '\t',stringsAsFactors = F)
  tmps<-tmps[,setdiff(1:ncol(tmps),c(grep('Asian',colnames(tmps)),grep('Eur',colnames(tmps))))]
  tmps_rnames<-rownames(tmps);  tmps1<-data.frame(tmps_rnames,tmps,stringsAsFactors = F)
  
  tmps2<-tmps1[as.character(tmps1$tmps_rnames)%in%genes_acmg,]
  tmps3_genes<-setdiff(genes_acmg,rownames(tmps2))
  tmps3<-data.frame(tmps3_genes,matrix(as.integer(0),ncol = ncol(tmps2)-1,nrow = length(tmps3_genes)),stringsAsFactors = F)
  colnames(tmps3)<-colnames(tmps2);  rownames(tmps3)<-tmps3_genes
  tmps4<-rbind(tmps3,tmps2)
  tmpplot2 <- tmps4
  
  tmps2<-tmps1[as.character(tmps1$tmps_rnames)%in%genes_top10p,]
  tmps3_genes<-setdiff(genes_top10p,rownames(tmps2))
  tmps3<-data.frame(tmps3_genes,matrix(as.integer(0),ncol = ncol(tmps2)-1,nrow = length(tmps3_genes)),stringsAsFactors = F)
  colnames(tmps3)<-colnames(tmps2);  rownames(tmps3)<-tmps3_genes
  tmps4<-rbind(tmps3,tmps2)
  tmps4 <- tmps4[setdiff(tmps4$tmps_rnames,tmpplot2$tmps_rnames),]
  tmpplot1 <- tmps4
  
  tmpplot <- rbind(tmpplot1,tmpplot2)
  
  tmptmp<-melt(tmpplot,id='tmps_rnames')
  tmptmp$tmps_rnames<-factor(tmptmp$tmps_rnames,levels = rev(c(genes_acmg,setdiff(genes_top10p,genes_acmg))))
  if(length(grep('AF',tmp[i]))!=0){tmptmp$value<-round(tmptmp$value,6)*100}#for if 
  p<-ggplot(tmptmp,aes(y=tmps_rnames,x=variable,fill=value))+geom_tile(color='grey',height=1)
  p<-p+geom_text(aes(label=ifelse(value!=0,value,NA),size=10))
  p<-p+theme_minimal()
  p<-p+theme(axis.text.x = element_text(size = 17,angle = 90,hjust = 1,vjust = 0.5),
             axis.text.y = element_text(size = 17),axis.title = element_blank(),legend.position = 'top',
             legend.text = element_text(size = 17,angle = 45,hjust = 1,vjust = 1),legend.title = element_blank(),
             panel.grid = element_blank())
  p<-p+scale_fill_gradient2(low = 'white',high = 'red')
  p
  ggsave(file.path(outdir,paste(tmp[i],'2.pdf',sep = '')),h=14,w=8)
  write.table(tmps1[genes_acmg,],file.path(outdir,paste(tmp[i],'2.txt',sep = '')),row.names = F,col.names = T,quote = F,sep = '\t')
}#for i
