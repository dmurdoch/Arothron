#' spherical.flipping
#'
#' Internal spherical flippping function
#' @param C numeric: coordinates of the point of view
#' @param mesh object of class mesh3d 
#' @param param1 numeric: first parameter for spherical flipping (usually ranged from 0.1 to 3, try!)
#' @param param2 numeric second paramter for spherical flipping (don't change it!)
#' @author Antonio Profico, Costantino Buzi, Marina Melchionna, Paolo Piras, Pasquale Raia, Alessio Veneziano
#' @references Profico, A., Schlager, S., Valoriani, V., Buzi, C., Melchionna, M., Veneziano, A., ... & Manzi, G. (2018). 
#' Reproducing the internal and external anatomy of fossil bones: Two new automatic digital tools. American Journal of Physical Anthropology, 166(4), 979-986.#' @export
#'
#' Katz, S., Tal, A., & Basri, R. (2007). Direct visibility of point sets. In ACM SIGGRAPH 2007 papers (pp. 24-es).
#' @export


spherical.flipping<-function(C,mesh,param1,param2){
C<-matrix(as.vector(C),ncol=3,nrow=1)
mesh<-mesh
P<-t(mesh$vb)[,-4]
numVer<-dim(P)[2]
numDim<-dim(P)[1]
P2<-P-(repmat(C,numDim,1))
normp<-rowSums(P2^2)
normp<-sqrt(normp)
param<-param1
R2<-matrix(repmat(max(normp)*(param2^param1),numDim, 1))
SF<-P2+2*repmat(R2-cbind(normp),1,3)*
P2/(repmat(cbind(normp),1,3))  
return(SF)
}
  
  