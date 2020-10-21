#' twodviews
#' Combine and calculate the PCscores matrix from a list of different landmark configurations to be combined
#' @param twodlist a list containing the landmark configurations of each anatomical view stored as separated lists
#' @param scale logical: TRUE for shape-space, FALSE for form-space 
#' @param vector numeric vector: defines which views are to be used
#' @return PCscores PC scores
#' @return PCs Pricipal Components (eigenvector matrix)
#' @return Variance table of the explained variance by the PCs
#' @return size vector containing the Centroid Size of each configuration
#' @return mshapes a list containing the mean shape of each landmark configuration
#' @return dims number of landmarks of each configuration
#' @return dimm dimension (2D or 3D) of each combined landmark configuration
#' @return twodlist the list used as input
#' @author Antonio Profico, Costantino Buzi, Marina Melchionna, Paolo Piras, Pasquale Raia, Alessio Veneziano
#' @references Profico, A., Piras, P., Buzi, C., Del Bove, A., Melchionna, M., Senczuk, G., ... & Manzi, G. (2019). 
#' Seeing the wood through the trees. Combining shape information from different landmark configurations. Hystrix, 157-165.
#' @examples
#' library(Morpho)
#' #load the 2D primate dataset
#' data("Lset2D_list")
#' length(Lset2D_list)
#' #combine the 2D datasets and PCA
#' combin2D<-twodviews(Lset2D_list,scale=TRUE,vector=c(1:5))
#' combin2D$size
#' #plot of the first two Principal Components
#' plot(combin2D$PCscores)
#' text(combin2D$PCscores,labels=rownames(combin2D$PCscores))
#' #load the 3D primate dataset
#' data("Lset3D_array")
#' #GPA and PCA
#' GPA_3D<-procSym(Lset3D_array)
#' #plot of the first two Principal Components
#' plot(GPA_3D$PCscores)
#' text(GPA_3D$PCscores,labels=rownames(GPA_3D$PCscores))
#' @export
#' 
twodviews<-function (twodlist, scale = TRUE, vector = NULL) {
  if(is.null(vector)){
    vector<-1:length(twodlist)
  }
  sizes = matrix(NA, ncol = length(vector), nrow = dim(listtoarray(twodlist[[1]]))[3])
  mshapes = list()
  dims = NULL
  Rots = NULL
  count <- 0
  dimm<-NULL
  for (h in vector) {
    count <- count + 1
    coo <- listtoarray(twodlist[[h]])
    dimm<-c(dimm,dim(coo)[2])
    dims <- c(dims, dim(coo)[1])
    if (scale == FALSE) {
      gpa <- procSym(coo, scale = FALSE, CSinit = FALSE)
    }
    else {
      gpa <- procSym(coo, scale = TRUE)
    }
    coo.r <- gpa$orpdata
    Rots <- cbind(Rots, (vecx(coo.r) * sqrt(dim(coo)[1] * dimm[h])))
    sizes[, count] <- gpa$size
    mshapes[[count]] <- gpa$mshape
  }
  
  Rots_pca <- prcomp(Rots, scale. = FALSE)
  values <- 0
  eigv <- Rots_pca$sdev^2
  values <- eigv[which(eigv > 1e-16)]
  lv <- length(values)
  PCs <- Rots_pca$rotation[, 1:lv]
  PCscores <- as.matrix(Rots_pca$x[, 1:lv])
  rownames(PCscores) <- names(twodlist[[1]])
  Variance <- cbind(sqrt(eigv), eigv/sum(eigv), cumsum(eigv)/sum(eigv)) * 100
  Variance <- Variance[1:lv, ]
  rotated <- list()
  for(i in 1:length(vector)){
    if(i ==1){
      rotated[[i]]<-vecx(Rots[,1:(dims[i]*dimm[i])],byrow = FALSE,revert = TRUE,lmdim = dimm[i])}else{
        rotated[[i]]<-vecx(Rots[,((dims[i-1]*dimm[i])+1):(sum(dims[(i-1):i])*dimm[i])],byrow = FALSE,revert = TRUE,lmdim = dimm[i])  
      }
  }
  
  out <- list(PCscores = PCscores, PCs = PCs, Variance = Variance, 
              size = sizes, rotated = rotated, mshapes = mshapes, dims = dims, 
              twodlist = twodlist, dimm=dimm)
  return(out)
}
