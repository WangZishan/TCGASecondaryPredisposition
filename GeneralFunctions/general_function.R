first_letter_up<-function(x){
  substr(x,1,1)<-toupper(substr(x,1,1))
  x
}#for function
#Transform a numeric vector into value of empirical cumulative distribution
myecdf<-function(x){
  if(sum(is.na(x))==length(x)){
    return(rep(NA,length(x)))
  }else{
    Fn<-ecdf(x)
    return(Fn(x))
  }
}#for myecdf
#round up function, which is slightly different with round function in R base package
round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}#function
#sum of not na element
sum_notna<-function(x){return(sum(x[!is.na(x)]))}#for function
sum_na<-function(x){return((sum(x[!is.na(x)])+sum(is.na(x))))}#for function
#return a logical vector showing whether each element is in duplicated element set
myduplicated<-function(x){return(x%in%x[duplicated(x)])}#for function

# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#     output_barplot
#
# output a simple file for barplot
#
# x - a numeric vector contain your interested information
# outfile - the name of file to save the plot
output_barplot<-function(x, stat = "count", title=NULL, xrange=NULL, outfile, h, w){
  data<-data.frame(x=x)
  p<-ggplot2::ggplot(data = data,aes(x=x))+ggplot2::geom_bar(stat = stat)
  p<-p+theme_bw()+
    theme(axis.text.x = element_text(size = 17),axis.text.y = element_text(size = 15),
          axis.title.x = element_blank(),axis.title.y = element_blank(),
          panel.grid.major = element_blank(),plot.margin = unit(c(0,0,0,0),"cm"),
          legend.text = element_text(size = 15,angle = 45),legend.title = element_text(size = 18),legend.position = "bottom")
  if(!is.null(xrange)){p<-p+xlim(xrange[1],xrange[2])}
  if(!is.null(title)){p<-p+labs(title = title)}
  p
  ggsave(outfile,h=h,w=w,useDingbat=F)
}#for function



