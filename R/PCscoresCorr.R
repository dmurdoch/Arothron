#' PCscoresCorr
#' Perform a correlation test between two matrices of PCscores
#' @param matrix1 matrix: first set of PC scores 
#' @param matrix2 matrix: second set of PC scores
#' @param nPCs numeric vector: specify which PC scores will be selected in the correlation test 
#' @return corr the correlation values associated to each pair of PC scores
#' @return p.values p-values associated to the correlation test
#' @author Antonio Profico, Costantino Buzi, Marina Melchionna, Paolo Piras, Pasquale Raia, Alessio Veneziano
#' @export
#' 
PCscoresCorr<-function(matrix1,matrix2,nPCs=1:5){
  Rcors<-NULL
  pvalues<-NULL
  for(i in nPCs){
    Rcorsi<-cor.test(matrix1[,i],matrix2[,i])$estimate
    if(Rcorsi<0){Rcors[i]<-cor.test(matrix1[,i],matrix2[,i]*-1)$estimate}else{
      Rcors[i]<-cor.test(matrix1[,i],matrix2[,i])$estimate
    }
    pvalues[i]<-cor.test(matrix1[,i],matrix2[,i])$p.value
  }  
  out<-list("corr"=Rcors,"p.values"=pvalues)
  return(out)
}
