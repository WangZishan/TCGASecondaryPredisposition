rm(list = ls()); tmp<-lapply(c('ggplot2','reshape2','grid'), library,character.only=T)
setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition/analysis/variant_distribution_of_ancestry/variants_pca_frq_nonCancer_pathogenic')

### load in required information of all kinds
# load in variant information
variant<-read.table('../../data/variant_data/variants_pca_frq_nonCancer_pathogenic.txt',header = T,sep = '\t',quote = '',stringsAsFactors = F)
# ancestry information
ancestry<-read.table('../../data/ancestry_info/ancestry_patient.txt',header = T,sep = '\t',stringsAsFactors = F)

### Pathogenic and likely pathogenic distribution across different ancestries
# count of spv among acmg classification across ancestries
tmp<-unique(variant[,c('HGVSg1','bcr_patient_barcode','ACMG','ancestry')])
p_count_long<-as.data.frame(table(tmp[,c('ACMG','ancestry')]),stringsAsFactors = F)
colnames(p_count_long)<-c('ACMG','Ancestry','Count'); p_count_long$label_loc<-NA
for (i in 1:nrow(p_count_long)) {p_count_long[i,'label_loc']<-sum(p_count_long[p_count_long$Ancestry==p_count_long[i,'Ancestry'],'Count'])}#for i
write.table(dcast(p_count_long,formula = Ancestry~ACMG,value.var = 'Count'),'Distribution_for_ACMG_classification-SPV_Count.txt',row.names = F,col.names = T,quote = F,sep = '\t')

p_count_long$Ancestry<-factor(p_count_long$Ancestry,levels = rev(c('Amr','Afr','Eur','Eas','Sas','Mix','Oth')))
p_count_long1<-p_count_long[p_count_long$ACMG=='likely pathogenic',]; p_count_long2<-p_count_long[p_count_long$ACMG=='pathogenic',]
plot_count<-ggplot()+geom_bar(data=p_count_long1,aes(x=Ancestry,y=label_loc),stat = 'identity',fill='#F8766D')
plot_count<-plot_count+geom_bar(data=p_count_long2,aes(x=Ancestry,y=Count),stat = 'identity',fill='#00BFC4')
plot_count<-plot_count+geom_text(data=p_count_long1,aes(x=Ancestry,y=label_loc*0.7,label=Count),hjust=0.5,size=5,colour='grey40')
plot_count<-plot_count+geom_text(data=p_count_long2,aes(x=Ancestry,y=Count/4,label=Count),hjust=0.5,size=5,colour='grey40')
plot_count<-plot_count+theme_bw()+  theme(axis.text = element_text(size = 20),axis.ticks.length.y = unit(0,'cm'),
                                          axis.title.x =  element_text(size = 23),axis.title.y = element_blank(),
                                          legend.title = element_text(size = 23),legend.text = element_text(size = 20),
                                          panel.grid.major.y = element_blank(),legend.position = 'none')
plot_count<-plot_count+labs()+coord_flip()+scale_y_continuous(trans = 'log10')+labs(y='SPV Count')
plot_count
ggsave('Distribution_for_ACMG_classification-SPV_Count.pdf',h=6,w=3,useDingbat=F)

#count/proportion of patient among acmg classification across ancestries
p_count_carrier_wide<-as.data.frame(matrix(nrow = 7,ncol = 3,dimnames = list(c(),c('Ancestry','likely pathogenic','pathogenic'))),stringsAsFactors = F)
p_count_carrier_wide$Ancestry<-c('Amr','Afr','Eur','Eas','Sas','Mix','Oth')
for (i in 1:nrow(p_count_carrier_wide)) {
  tmpvariant<-variant[variant[,'ancestry']==p_count_carrier_wide[i,'Ancestry'],]
  p_count_carrier_wide[i,'pathogenic']<-length(unique(tmpvariant[tmpvariant[,'ACMG']=='pathogenic','bcr_patient_barcode']))
  p_count_carrier_wide[i,'likely pathogenic']<-length(unique(setdiff(tmpvariant[tmpvariant[,'ACMG']=='likely pathogenic','bcr_patient_barcode'],tmpvariant[tmpvariant[,'ACMG']=='pathogenic','bcr_patient_barcode'])))
}#for i
tmp<-p_count_carrier_wide[,2:3]; tmp<-tmp/as.numeric(table(ancestry$ancestry)[as.character(p_count_carrier_wide$Ancestry)])
p_percentage_carrier_wide<-data.frame(p_count_carrier_wide$Ancestry,tmp,stringsAsFactors = F);colnames(p_percentage_carrier_wide)<-c('Ancestry',gsub('\\.',' ',colnames(p_percentage_carrier_wide)[2:3]))
write.table(cbind(p_percentage_carrier_wide,p_count_carrier_wide),'Distribution_for_ACMG_classification-SPV_carrier_percentage.txt',row.names = F,col.names = T,quote = F,sep = '\t')

p_percentage_carrier_long<-melt(p_percentage_carrier_wide,id.vars = 'Ancestry'); colnames(p_percentage_carrier_long)<-c('Ancestry','ACMG','Percentage')
p_percentage_carrier_long$Percentage<-round(p_percentage_carrier_long$Percentage*100,2)
p_percentage_carrier_long$label_loc<-NA
for (i in 1:nrow(p_percentage_carrier_long)) {
  tmpp_percentage_carrier_long<-p_percentage_carrier_long[p_percentage_carrier_long$Ancestry==p_percentage_carrier_long[i,'Ancestry'],]
  p_percentage_carrier_long[i,'label_loc']<-sum(tmpp_percentage_carrier_long$Percentage)
}#for i

p_percentage_carrier_long$Ancestry<-factor(p_percentage_carrier_long$Ancestry,levels = rev(c('Amr','Afr','Eur','Eas','Sas','Mix','Oth')))
p_percentage_carrier_long1<-p_percentage_carrier_long[p_percentage_carrier_long$ACMG=='likely pathogenic',]
p_percentage_carrier_long2<-p_percentage_carrier_long[p_percentage_carrier_long$ACMG=='pathogenic',]
plot_percentage<-ggplot()+geom_bar(data=p_percentage_carrier_long1,aes(x=Ancestry,y=label_loc),stat = 'identity',fill='#F8766D')
plot_percentage<-plot_percentage+geom_bar(data = p_percentage_carrier_long2,aes(x=Ancestry,y=Percentage),stat = 'identity',fill='#00BFC4')
plot_percentage<-plot_percentage+geom_text(data=p_percentage_carrier_long1,aes(x=Ancestry,y=label_loc*0.9,label=Percentage),hjust=0.5,size=5,colour='grey40')
plot_percentage<-plot_percentage+geom_text(data = p_percentage_carrier_long2,aes(x=Ancestry,y=Percentage/2,label=Percentage),hjust=0.5,size=5,colour='grey40')
plot_percentage<-plot_percentage+theme_bw()+theme(axis.text = element_text(size = 20),axis.ticks.length.y = unit(0,'cm'),
                                                  axis.title.x = element_text(size = 23),axis.title.y = element_blank(),
                                                  legend.title = element_text(size = 23),legend.text = element_text(size = 20),
                                                  panel.grid.major.y = element_blank(),legend.position = 'none')
plot_percentage<-plot_percentage+scale_fill_discrete(name='Variant Classification')
plot_percentage<-plot_percentage+labs(y='Percentage (%)')+coord_flip()+scale_y_reverse()
plot_percentage
ggsave('Distribution_for_ACMG_classification-SPV_carrier_percentage.pdf',h=6,w=3,useDingbat=F)
