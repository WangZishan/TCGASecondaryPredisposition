column_average<-function(x,y){
  ydup<-unique(y[duplicated(y)])
  ydup_l<-length(ydup)
  
  tmpind<-y%in%ydup
  xuni<-x[,!tmpind]
  xdup<-x[,tmpind]
  xdup_colnames<-y[tmpind]
  
  xdup_uni<-matrix(nrow = nrow(x),ncol = ydup_l)
  colnames(xdup_uni)<-ydup
  rownames(xdup_uni)<-rownames(x)
  for (i in 1:ydup_l) {
    tmpxdup_uni<-xdup[,xdup_colnames==ydup[i]]
    xdup_uni[,i]<-rowMeans(tmpxdup_uni)
  }#for i
  
  x_merge<-data.frame(xdup_uni,xuni,stringsAsFactors = F)
  return(x_merge)
}#for function
