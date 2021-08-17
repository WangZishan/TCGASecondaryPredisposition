rm(list = ls()); tmp<-lapply(c('ggplot2','scatterpie'), library, character.only=T)
setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition')

### load in related info
# sample-variant info
variant<-read.table('analysis/data/variant_data/variants_pca_frq_nonCancer_pathogenic.txt',header = T,sep = '\t',stringsAsFactors = F,quote = '')
# ASE analysis result info
ase<-read.table('analysis/allele_expression/SPV/ASE_basic-allele_count.txt',header = T,sep = '\t',stringsAsFactors = F)
# genes of interest
tmpgenes<-read.table('analysis/allele_expression/SPV/ASE_Gene_Enrichment.txt',header = T,sep = '\t',stringsAsFactors = F)
genes<-c(tmpgenes$Symbol_Merge[tmpgenes$Classification=='Significant'],tmpgenes$Symbol_Merge[tmpgenes$Classification=='Suggestive'])

### Extract sample-variant info of genes of interest
# merge sample-variant info and ASE analysis info
tmp<-ase[match(paste(variant$HGVSg1,variant$bcr_patient_barcode),paste(ase$HGVSg1,ase$bcr_patient_barcode)),]
colnames(tmp)<-paste('ASE',colnames(ase),sep = '_')
variant<-cbind(variant,tmp)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# genes of intersts must have no overlap with variant corresponding to multiple gene, otherwise the code need to modify a little
#------------------------------------------------------------------------------------------------------------------------------------------------------
table(unique(unlist(strsplit(variant$Symbol_Merge[grep(';',variant$Symbol_Merge)],';')))%in%genes)

# extract process
#------------------------------------------------------------------------------------------------------------------------------------------------------
# manually construct the typ_ma matrix
# only start and stop position of the following variants should be equal
#------------------------------------------------------------------------------------------------------------------------------------------------------
variant<-variant[variant$Symbol_Merge%in%genes,]; unique(variant$Function); table(variant$Start==variant$Stop)
typ_ma<-data.frame(c('stop-gain','missense','splice-acceptor','synonymous'),c('N','M','L','S'),stringsAsFactors = F)
for (i in 1:length(genes)) {
  tmpvariant<-variant[variant$Symbol_Merge==genes[i],]; typ<-typ_ma[match(tmpvariant$Function,typ_ma[,1]),2]
  tmpvariant1<-data.frame(paste(tmpvariant$HGVSg1,tmpvariant$HGVSp,sep = '_'),paste('chr',tmpvariant$Chromosome,':',tmpvariant$Start,sep = ''),typ,stringsAsFactors = F)
  tmpvariant2<-data.frame(tmpvariant$HGVSg1,paste('chr',tmpvariant$Chromosome,':',tmpvariant$Start,sep = ''),typ,tmpvariant[,c('Function','Feature','ASE_pval_fdr','ASE_class','HGVSp')],stringsAsFactors = F)
  write.table(tmpvariant1,file.path('analysis/allele_expression/SPV/',paste('GeneVariant_Distribution-',genes[i],'.txt',sep = '')),row.names = F,col.names = F,quote = F,sep = ';')
  write.table(tmpvariant2,file.path('analysis/allele_expression/SPV/',paste('GeneVariant_Distribution-',genes[i],'_transcript.txt',sep = '')),row.names = F,col.names = F,quote = F,sep = ';')
}#for i

### Scatterpie plot
m <- matrix(nrow = 0,ncol = 10,dimnames = list(NULL,c('Gene','GeneIndex','HGVSp','HGVSpIndex','Index','Count','Sig','Sug','NotASE','Insufficient')))
for (i in 1:length(genes)) {
  tmpvariant<-variant[variant$Symbol_Merge==genes[i],]
  hgvsp<-unique(tmpvariant$HGVSp)
  for (j in 1:length(hgvsp)) {
    tmptmpvariant <- tmpvariant[tmpvariant$HGVSp==hgvsp[j],]
    tmpm <- matrix(nrow = 1,ncol = 10,dimnames = list(NULL,c('Gene','GeneIndex','HGVSp','HGVSpIndex','Index','Count','Sig','Sug','NotASE','Insufficient')))
    tmpm[1,'Gene'] <- genes[i]; tmpm[1,"GeneIndex"] <- i
    tmpm[1,"HGVSp"] <- hgvsp[j]; tmpm[1,"HGVSpIndex"] <- j
    tmpm[1,"Count"] <- nrow(tmptmpvariant)
    tmpm[1,"Sig"] <- sum(tmptmpvariant$ASE_class=='Significant',na.rm = T)
    tmpm[1,"Sug"] <- sum(tmptmpvariant$ASE_class=='Suggestive',na.rm = T)
    tmpm[1,"NotASE"] <- sum(tmptmpvariant$ASE_class=='Not ASE',na.rm = T)
    tmpm[1,"Insufficient"] <- sum(tmptmpvariant$ASE_class=='Insufficient read count data',na.rm = T)+sum(is.na(tmptmpvariant$ASE_class))
    m <- rbind(m,tmpm)
  }# for j
}# for i
m[,"Index"] <- 1:nrow(m)
m<-as.data.frame(m,stringsAsFactors = F)
for (i in c('GeneIndex','HGVSpIndex','Index','Count','Sig','Sug','NotASE','Insufficient')) {
  m[,i] <- as.numeric(m[,i])
}# for
p <- ggplot() + geom_scatterpie(aes(x=GeneIndex, y=HGVSpIndex, group=Index, r=0.1+Count*0.02), data=m,
                                cols=c('Sig','Sug','NotASE','Insufficient'), color=NA) + coord_equal()
p <- p + scale_fill_manual(values = c('Sig'='#6159A4','Sug'='#CE4A08','NotASE'='#1D8F64','Insufficient'='#B3B3B3'))
p <- p + theme_bw()
p
ggsave(p,filename = 'analysis/allele_expression/SPV/GeneVariant_Distribution-Pie.pdf',h=10,w=10,useDingbat=F)
writexl::write_xlsx(m,'analysis/allele_expression/SPV/GeneVariant_Distribution-Pie.xlsx',format_headers = F)
