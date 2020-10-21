#' CScorreffect
#' Plot showing the correlation in the shape space between original and combined dataset omitting or including the normalization factors calculated with Arothron and MLECScorrection 
#' @param array1 array: first set of landmark configuration
#' @param array2 array: second set of landmark configuration
#' @param nPCs numeric vector: specify which PC scores will be selected in the correlation test 
#' @param from numeric: the lower interval of the normalization factor distribution
#' @param to numeric: the lower interval of the normalization factor distribution
#' @param length.out numeric: number of values ranged between from and to 
#' @return PCscores PCscores matrix of the combined dataset applying the normalization factor calculated by using the maximum likelihood estimation
#' @return PCs PCs matrix of the combined dataset applying the normalization factor calculated by using the maximum likelihood estimation
#' @return corr mean correlation between original and combined dataset
#' @return CSratios normalization factor calculated by using the maximum likelihood estimation
#' @author Antonio Profico, Costantino Buzi, Marina Melchionna, Paolo Piras, Pasquale Raia, Alessio Veneziano
#' @examples
#' \dontrun{
#' # Femora case study
#' data(femsets)
#' all_pois<-matrix(1:(200*61),nrow=61,ncol=200,byrow = FALSE)
#' set_ext_100<-femsets[all_pois[,1:100],,]
#' 
#' set_int_100<-femsets[all_pois[,101:200],,]
#' set_int_50<-set_int_100[c(matrix(1:6100,ncol=61)[seq(1,100,2),]),,]
#' set_int_20<-set_int_100[c(matrix(1:6100,ncol=61)[seq(1,100,5),]),,]
#' set.seed(123)
#' sel<-sample(1:100,10)
#' set_int_10r<-set_int_100[c(matrix(1:6100,ncol=61)[sel,]),,]
#' 
#' CScorreffect(set_ext_100,set_int_50,nPCs=1:3)
#' CScorreffect(set_ext_100,set_int_20,nPCs=1:3)
#' CScorreffect(set_ext_100,set_int_10r,nPCs=1:3)
#' }
#' @export
 




CScorreffect<-function(array1,array2,nPCs=c(1:3),from=0.02,to=0.90,length.out= 100){
  
  data1<-MLECScorrection(array1,array2)
  Array12<-bindArr(array1,array2,along=1)
  PCscores_ref<-procSym(Array12,scale=TRUE)$PCscores
  seqs<-seq(from=from,to=to,length.out = length.out)
  gpa1 <- procSym(array1, scale = TRUE)
  gpa2 <- procSym(array2, scale = TRUE)
  CSratios<-NULL
  corj<-NULL
  for(j in 1:length(seqs)){
    coo.r1<-gpa1$orpdata*seqs[j]
    coo.r2<-gpa2$orpdata*(1-seqs[j])
    Rots <- cbind(vecx(coo.r1), vecx(coo.r2))
    Rots_pca <- prcomp(Rots, scale. = FALSE)
    values <- 0
    eigv <- Rots_pca$sdev^2
    values <- eigv[which(eigv > 1e-16)]
    lv <- length(values)
    PCs <- Rots_pca$rotation[, 1:lv]
    PCscores <- as.matrix(Rots_pca$x[, 1:lv])
    corj[j]<-mean(PCscoresCorr(PCscores_ref,PCscores,nPCs)$corr)
    CSratios[j]<-(seqs[j]/(1-seqs[j]))
  }
  plot(CSratios,corj,type="b",pch=19,main="red=Arothron CS correction, black= MLE CS correction,blue=no CS correction",
       xlab="Centroid size correction", ylab="Correlation")
  abline(h=data1$corr,lty=2,lwd=2)
  abline(v=data1$CScorr,lty=2,lwd=2)
  datalist1<-list("ext"=arraytolist(array1),"int"=arraytolist(array2))
  aro_meth<-twodviews(datalist1,scale=TRUE,vector = 1:2)
  aro_res<-mean(PCscoresCorr(PCscores_ref,aro_meth$PCscores,nPCs = nPCs)$corr)
  abline(h=aro_res,lwd=2,lty=3,col="red")
  aro_CSratio<- sqrt(dim(array1)[1] * 3)/sqrt(dim(array2)[1] * 3)
  abline(v=aro_CSratio,lwd=2,lty=3,col="red")
  
  coo.r1<-gpa1$orpdata
  coo.r2<-gpa2$orpdata
  Rots <- cbind(vecx(coo.r1), vecx(coo.r2))
  Rots_pca <- prcomp(Rots, scale. = FALSE)
  values <- 0
  eigv <- Rots_pca$sdev^2
  values <- eigv[which(eigv > 1e-16)]
  lv <- length(values)
  PCs <- Rots_pca$rotation[, 1:lv]
  PCscores <- as.matrix(Rots_pca$x[, 1:lv])
  
  standcorr<-mean(PCscoresCorr(PCscores_ref,PCscores,nPCs)$corr)
  noCSratio<- mean(apply(array1,3,cSize))/mean(apply(array2,3,cSize))
  abline(v=1,lwd=2,lty=2,col="blue")
  abline(h=standcorr,lwd=2,lty=2,col="blue")
  
  out<-list("MLE-cor"=data1$corr,"MLE-CS"=data1$CScorr,
            "Aro-cor"=aro_res,"Aro-CS"=aro_CSratio,
            "Sta-cor"=standcorr,"Sta-CS"=1)
  return(out)
}
