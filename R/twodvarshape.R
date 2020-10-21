#' twodvarshape 
#' Calculates the shape variation associated to a value of PC scores associated to a specific combined landmark configuration or view
#' @param twodviews_ob object from twodviews()
#' @param scores numeric: the values of the PC scores for which the visualization is called
#' @param PC PC chosen
#' @param view numeric: which landmark configuration will be used to build the shape variation 
#' @return mat matrix of coordinates associated to the called shape variation
#' @author Antonio Profico, Costantino Buzi, Marina Melchionna, Paolo Piras, Pasquale Raia, Alessio Veneziano
#' @references Profico, A., Piras, P., Buzi, C., Del Bove, A., Melchionna, M., Senczuk, G., ... & Manzi, G. (2019). 
#' Seeing the wood through the trees. Combining shape information from different landmark configurations. Hystrix, 157-165.
#' @examples
#' library(Arothron)
#' #load the 2D primate dataset
#' data("Lset2D_list")
#' #combine the 2D datasets and PCA
#' combin2D<-twodviews(Lset2D_list,scale=TRUE,vector=c(1:5))
#' #calculate the shape variation associated to the negative extreme value of PC1
#' min_PC1<-twodvarshape(combin2D,min(combin2D$PCscores[,1]),1,5)
#' plot(min_PC1,asp=1)
#' #calculate the shape variation associated to the positive extreme value of PC1
#' max_PC1<-twodvarshape(combin2D,max(combin2D$PCscores[,1]),1,5)
#' plot(max_PC1,asp=1)
#' @export

twodvarshape<-function (twodviews_ob, scores, PC, view) 
{
  pos_pcs <- twodviews_ob$dims * twodviews_ob$dimm
  if (view == 1) {
    sel_pcs <- 1:pos_pcs[view]
  }
  if (view == (length(twodviews_ob$dims))) {
    sel_pcs <- (sum(pos_pcs[1:(view - 1)]) + 1):sum(pos_pcs[1:view])
  }
  if (view != 1 & view != (length(view))) {
    sel_pcs <- (sum(pos_pcs[1:(view - 1)]) + 1):(sum(pos_pcs[1:(view - 
                                                                  1)]) + pos_pcs[view])
  }
  mshape <- twodviews_ob$mshapes[[view]] * sqrt(twodviews_ob$dims[view] * twodviews_ob$dimm[view])
  PCs <- twodviews_ob$PCs[sel_pcs, PC]
  mat <- showPC(scores, PCs, mshape)
  return(mat)
}

