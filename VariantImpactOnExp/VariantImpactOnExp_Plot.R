rm(list = ls()); tmp<-lapply(c('ggplot2','ggrepel'),library,character.only=T)
setwd('C:/Users/wangz21/Box Sync/TCGASecondaryPredisposition')
source('analysis/functions/global_aes_out.R')
Enriched_status<-c('Sig.'='#F8766D','Sug.'='#00BA38','Not Sig.'='#619CFF')
outdir<-indir<-file.path(getwd(),'analysis','VariantImpactOnExp','exp_percentile')

### Load in related info
# linear regression model result
individual_cancer<-read.table(file.path(indir,'VariantImpactOnExp-individual_cancer.txt'),header = T,sep = '\t',stringsAsFactors = F)
individual_normal<-read.table(file.path(indir,'VariantImpactOnExp-individual_normal.txt'),header = T,sep = '\t',stringsAsFactors = F)
pan_cancer<-read.table(file.path(indir,'VariantImpactOnExp-pan_cancer.txt'),header = T,sep = '\t',stringsAsFactors = F)
pan_normal<-read.table(file.path(indir,'VariantImpactOnExp-pan_normal.txt'),header = T,sep = '\t',stringsAsFactors = F)

### Statistics of significant linear regression model result
m<-matrix(nrow = 4,ncol = 6,dimnames = list(c('individual_cancer','individual_normal','pan_cancer','pan_normal'),c('FDR_0.15','FDR_0.05','FDR_0.01','pvalue_0.15','pvalue_0.05','pvalue_0.01')))
tmp<-individual_cancer$FDR1; tmp<-tmp[!is.na(tmp)]
m[1,3]<-sum(tmp<0.01); m[1,2]<-sum(tmp<0.05); m[1,1]<-sum(tmp<0.15)
tmp<-individual_normal$FDR1; tmp<-tmp[!is.na(tmp)]
m[2,3]<-sum(tmp<0.01); m[2,2]<-sum(tmp<0.05); m[2,1]<-sum(tmp<0.15)
tmp<-pan_cancer$FDR1; tmp<-tmp[!is.na(tmp)]
m[3,3]<-sum(tmp<0.01); m[3,2]<-sum(tmp<0.05); m[3,1]<-sum(tmp<0.15)
tmp<-pan_normal$FDR1; tmp<-tmp[!is.na(tmp)]
m[4,3]<-sum(tmp<0.01); m[4,2]<-sum(tmp<0.05); m[4,1]<-sum(tmp<0.15)
tmp<-individual_cancer$P_Chi1; tmp<-tmp[!is.na(tmp)]
m[1,6]<-sum(tmp<0.01); m[1,5]<-sum(tmp<0.05); m[1,4]<-sum(tmp<0.15)
tmp<-individual_normal$P_Chi1; tmp<-tmp[!is.na(tmp)]
m[2,6]<-sum(tmp<0.01); m[2,5]<-sum(tmp<0.05); m[2,4]<-sum(tmp<0.15)
tmp<-pan_cancer$P_Chi1; tmp<-tmp[!is.na(tmp)]
m[3,6]<-sum(tmp<0.01); m[3,5]<-sum(tmp<0.05); m[3,4]<-sum(tmp<0.15)
tmp<-pan_normal$P_Chi1; tmp<-tmp[!is.na(tmp)]
m[4,6]<-sum(tmp<0.01); m[4,5]<-sum(tmp<0.05); m[4,4]<-sum(tmp<0.15)
write.table(m,file.path(outdir,'VariantImpactOnExp_Plot-Sta.txt'),row.names = T,col.names = T,quote = F,sep = '\t')

### Volcano plot of significant genes in the above linear regression model result
individual_cancer1<-individual_cancer[!is.na(individual_cancer$P_Chi1),]
if(nrow(individual_cancer1)!=0){
  p <- ggplot(data=individual_cancer1,aes(y=coefficient,x=cancer,color = cancer))
  p <- p + geom_point(aes(size=-log10(P_Chi1),stroke=0),alpha=0.8)
  p <- p + getPCACancerColor()
  p <- p + geom_text(aes(label=ifelse(P_Chi1<0.15&P_Chi1>=0.05,gene,NA),size=1),color='black',nudge_y = -0.015,show.legend = F)
  p <- p + geom_text(aes(label=ifelse(P_Chi1<0.05&P_Chi1>=0.01,gene,NA),size=1),color='blue',nudge_y = -0.015,show.legend = F)
  p <- p + geom_text(aes(label=ifelse(P_Chi1<0.01,gene,NA),size=1),color='red',nudge_y = -0.015,show.legend = F)
  p <- p + labs(x='cancer',y='coefficient')
  p <- p + geom_hline(yintercept = 0, alpha=0.5)
  p <- p + theme_bw() + 
    theme(axis.text.x = element_text(colour = 'black', size = 8, angle = 90, vjust = 0.5),axis.text.y = element_text(colour = 'black',size = 8),axis.ticks = element_blank())
  p <- p + guides(fill = guide_legend(override.aes = aes(label='')))
  p
  ggsave(file.path(outdir,'VariantImpactOnExp_Plot-individual_cancer-pvalue_0.05.pdf'),h=7,w=8,useDingbat=F)
  
  p<-ggplot(data = individual_cancer1,aes(x=coefficient,y=-log10(P_Chi1),color=cancer))+
    geom_point(alpha=0.8,aes(size=num_mutation_patient/num_patient,stroke=0))
  p <- p + geom_label_repel(aes(label=ifelse(P_Chi1<0.15,paste(cancer,gene,sep = '\n'),NA)),show.legend =FALSE,size=5)
  p <- p + getPCACancerColor()
  p<-p+geom_vline(xintercept = 0,linetype='dashed')+geom_hline(yintercept = -log10(0.05),linetype='dashed')+geom_hline(yintercept = -log10(0.15),linetype='dashed')
  p<-p+theme_bw()+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 18),
          legend.title = element_text(size = 20),legend.text = element_text(size = 18))
  p<-p+labs(size='Carrier\nproportion',x='Coefficient',y='-log10(P value)')
  p
  ggsave(file.path(outdir,'VariantImpactOnExp_Plot-individual_cancer-pvalue_0.05_volcano.pdf'),h=7,w=8,useDingbat=F)
}

individual_normal1<-individual_normal[!is.na(individual_normal$P_Chi1),]
if(nrow(individual_normal1)!=0){
  p <- ggplot(data=individual_normal1,aes(y=coefficient,x=cancer,color = cancer))
  p <- p + geom_point(aes(size=-log10(P_Chi1),stroke=0),alpha=0.8)
  p <- p + getPCACancerColor()
  p <- p + geom_text(aes(label=ifelse(P_Chi1<0.15&P_Chi1>=0.05,gene,NA),size=1),color='black',nudge_y = -0.015,show.legend = F)
  p <- p + geom_text(aes(label=ifelse(P_Chi1<0.05&P_Chi1>=0.01,gene,NA),size=1),color='blue',nudge_y = -0.015,show.legend = F)
  p <- p + geom_text(aes(label=ifelse(P_Chi1<0.01,gene,NA),size=1),color='red',nudge_y = -0.015,show.legend = F)
  p <- p + labs(x='cancer',y='coefficient')
  p <- p + geom_hline(yintercept = 0, alpha=0.5)
  p <- p + theme_bw() + 
    theme(axis.text.x = element_text(colour = 'black', size = 8, angle = 90, vjust = 0.5),axis.text.y = element_text(colour = 'black',size = 8),axis.ticks = element_blank())
  p <- p + guides(fill = guide_legend(override.aes = aes(label='')))
  p
  ggsave(file.path(outdir,'VariantImpactOnExp_Plot-individual_normal-pvalue_0.05.pdf'),h=7,w=8,useDingbat=F)
  
  p<-ggplot(data = individual_normal1,aes(x=coefficient,y=-log10(P_Chi1),color=cancer))+
    geom_point(alpha=0.8,aes(size=num_mutation_patient/num_patient,stroke=0))
  p <- p + geom_label_repel(aes(label=ifelse(P_Chi1<0.15,paste(cancer,gene,sep = '\n'),NA)),show.legend =FALSE,size=5)
  p <- p + getPCACancerColor()
  p<-p+geom_vline(xintercept = 0,linetype='dashed')+geom_hline(yintercept = -log10(0.05),linetype='dashed')+geom_hline(yintercept = -log10(0.15),linetype='dashed')
  p<-p+theme_bw()+
    theme(axis.title = element_text(size = 20),axis.text = element_text(size = 18),
          legend.title = element_text(size = 20),legend.text = element_text(size = 18))
  p<-p+labs(size='Carrier\nproportion',x='Coefficient',y='-log10(P value)')
  p
  ggsave(file.path(outdir,'VariantImpactOnExp_Plot-individual_normal-pvalue_0.05_volcano.pdf'),h=7,w=8,useDingbat=F)
}

pan_normal1<-pan_normal[!is.na(pan_normal$FDR1),]
tmp<-rep('Not Sig.',nrow(pan_normal1))
tmp[pan_normal1$FDR1<0.15 & pan_normal1$FDR1>=0.05]<-'Sug.'
tmp[pan_normal1$FDR1<0.05]<-'Sig.'
pan_normal1$FDR1_class<-tmp
if(nrow(pan_normal1)!=0){
  p<-ggplot(data=pan_normal1,aes(x=coefficient,y=-log10(FDR1)))+
    geom_point(alpha=0.8,aes(size=num_mutation_patient/num_patient,stroke=0))
  p<-p+scale_color_manual(values = Enriched_status)
  p<-p+geom_label_repel(aes(label=ifelse(FDR1_class!='Not Sig.',gene,NA)),show.legend =FALSE,size=5)
  p<-p+theme_bw()
  p<-p+theme(axis.title = element_text(size = 20),axis.text =  element_text(size = 17),
             legend.title = element_text(size = 20),legend.text = element_text(size = 17),
             panel.grid.major = element_line(color = 'grey95'),panel.grid.minor = element_line(color = 'grey95'),
             plot.margin = unit(c(0,0,0,0),'cm'))
  p<-p+geom_vline(xintercept = 0,col='white')+geom_hline(yintercept = -log10(0.05),col='white')+geom_hline(yintercept = -log10(0.15),col='white')
  p<-p+geom_vline(xintercept = 0,linetype='dashed')+geom_hline(yintercept = -log10(0.05),linetype='dashed')+geom_hline(yintercept = -log10(0.15),linetype='dashed')
  p
  ggsave(file.path(outdir,'VariantImpactOnExp_Plot-pan_normal-FDR_0.05_volcano.pdf'),h=7, w=8,useDingbat=F)
}

pan_cancer1<-pan_cancer[!is.na(pan_cancer$FDR1),]
tmp<-rep('Not Sig.',nrow(pan_cancer1))
tmp[pan_cancer1$FDR1<0.15 & pan_cancer1$FDR1>=0.05]<-'Sug.'
tmp[pan_cancer1$FDR1<0.05]<-'Sig.'
pan_cancer1$FDR1_class<-tmp
if(nrow(pan_cancer1)!=0){
  p<-ggplot(data=pan_cancer1,aes(x=coefficient,y=-log10(FDR1),color=FDR1_class))+
    geom_point(alpha=0.8,aes(size=num_mutation_patient/num_patient,stroke=0))
  p<-p+scale_color_manual(values = Enriched_status)
  p<-p+geom_label_repel(aes(label=ifelse(FDR1_class!='Not Sig.',gene,NA)),show.legend =FALSE,size=5)
  p<-p+theme_bw()
  p<-p+theme(axis.title = element_text(size = 20),axis.text =  element_text(size = 17),axis.text.x = element_text(hjust = 0.8),
             legend.title = element_text(size = 20),legend.text = element_text(size = 17),
             panel.grid.major = element_line(color = 'grey95'),panel.grid.minor = element_line(color = 'grey95'),
             plot.margin = unit(c(0,0,0,0),'cm'))
  p<-p+geom_vline(xintercept = 0,col='white')+geom_hline(yintercept = -log10(0.05),col='white')+geom_hline(yintercept = -log10(0.15),col='white')
  p<-p+geom_vline(xintercept = 0,linetype='dashed')+geom_hline(yintercept = -log10(0.05),linetype='dashed')+geom_hline(yintercept = -log10(0.15),linetype='dashed')
  p<-p+labs(x='Coefficient',y='-log10(FDR)',color='Sig. or not Sig.',size='Carrier\nproportion')
  p
  ggsave(file.path(outdir,'VariantImpactOnExp_Plot-pan_cancer-FDR_0.05_volcano.pdf'),h=7, w=8,useDingbat=F)
}
