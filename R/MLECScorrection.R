#' MLECScorrection
#' Maximum Likelihood Estimation of the normalization factor to be applied to optimize the correlation between two landmark configurations to be combined by using twodviews and 
#' @param array1 array: first set of landmark configuration
#' @param array2 array: second set of landmark configuration
#' @param scale logical: if FALSE the analysis is performed in the shape space, if TRUE the analysis is performed in the size and shape space (gpa without scaling) 
#' @param nPCs numeric vector: specify which PC scores will be selected in the correlation test 
#' @return PCscores PCscores matrix of the combined dataset applying the normalization factor calculated by using the maximum likelihood estimation
#' @return PCs PCs matrix of the combined dataset applying the normalization factor calculated by using the maximum likelihood estimation
#' @return corr mean correlation between original and combined dataset
#' @return CSratios normalization factor calculated by using the maximum likelihood estimation
#' @author Antonio Profico, Costantino Buzi, Marina Melchionna, Paolo Piras, Pasquale Raia, Alessio Veneziano
#' @export
#' 

MLECScorrection<-function(array1,array2,scale=TRUE,nPCs=1:5){ 

  Array12<-bindArr(array1,array2,along=1)
  if (scale == TRUE) {
    PCscores_ref<-procSym(Array12,scale=TRUE)$PCscores}else{
      PCscores_ref<-procSym(Array12,scale=FALSE,CSinit = FALSE)$PCscores
    }
  
  
  if (scale == FALSE) {
    gpa1 <- procSym(array1, scale = FALSE, CSinit = FALSE)
    gpa2 <- procSym(array2, scale = FALSE, CSinit = FALSE)}else {
      gpa1 <- procSym(array1, scale = TRUE)
      gpa2 <- procSym(array2, scale = TRUE)
    }
  
  foo<-function(ssq){
    coo.r1<-gpa1$orpdata* ssq
    coo.r2<-gpa2$orpdata*(1-ssq)
    Rots <- cbind(vecx(coo.r1), vecx(coo.r2))
    Rots_pca <- prcomp(Rots, scale. = FALSE)
    values <- 0
    eigv <- Rots_pca$sdev^2
    values <- eigv[which(eigv > 1e-16)]
    lv <- length(values)
    PCs <- Rots_pca$rotation[, 1:lv]
    PCscores <- as.matrix(Rots_pca$x[, 1:lv])
    corj<-mean(PCscoresCorr(PCscores_ref,PCscores,nPCs)$corr)   
    abs(1-corj)
  }
  
  
  # fit<-mle(foo, start = list(ssq=0.5), lower=.001, upper=0.999, method = "L-BFGS-B" )
  fit<-mle(foo, start = list(ssq=0.1), lower=.001, upper=0.999, method = "L-BFGS-B" )
  
  fit@coef->SSQ
  
  
  
  coo.r1<-gpa1$orpdata*SSQ
  coo.r2<-gpa2$orpdata*(1-SSQ)
  Rots <- cbind(vecx(coo.r1), vecx(coo.r2))
  Rots_pca <- prcomp(Rots, scale. = FALSE)
  values <- 0
  eigv <- Rots_pca$sdev^2
  values <- eigv[which(eigv > 1e-16)]
  lv <- length(values)
  PCs <- Rots_pca$rotation[, 1:lv]
  PCscores <- as.matrix(Rots_pca$x[, 1:lv])
  corj<-mean(PCscoresCorr(PCscores_ref,PCscores,nPCs)$corr)
  CSratios<-(SSQ/(1-SSQ))
  
  out<-list("PCscores"=PCscores,"PCs"=PCs,"corr"=corj,"CScorr"=CSratios)
  return(out)
}
