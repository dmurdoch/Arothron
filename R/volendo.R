#' volendo
#'
#' Calculate the volume of a mesh by using a voxel-based method
#' @param mesh object of class mesh3d
#' @param alpha_vol numeric: alpha shape for construction external concave hull
#' @param ncells numeric: approximative number of cell for 3D grid construction
#' @return vol numeric: volume of the mesh expressed in cc
#' @author Antonio Profico, Costantino Buzi, Marina Melchionna, Paolo Piras, Pasquale Raia, Alessio Veneziano
#' @examples
#' \dontrun{
#' #load the human skull
#' library(rgl)
#' data(human_skull)
#' sapendo<-endomaker(human_skull,param1_endo = 1.0,vol=FALSE, num.cores=NULL)
#' volsap<-volendo(sapendo$endocast)
#' }
#' @export

volendo<-function (mesh, alpha_vol = 100, ncells = 1e+05) 
{
   if (checkFaceOrientation(mesh) == TRUE){
    mesh<-invertFaces(mesh)
  } 
  
  conv_endo <- ashape3d(vert2points(mesh), alpha = alpha_vol, 
                        pert = TRUE, eps = 1e-06)
  x <- conv_endo
  selrows = x$triang[, 8 + 1] %in% 2L
  tr <- x$triang[selrows, c("tr1", "tr2", "tr3")]
  m = rgl::tmesh3d(vertices = t(x$x), indices = t(tr), homogeneous = FALSE)
  sn = alphashape3d::surfaceNormals(x, indexAlpha = 1)
  fn = facenormals(m)
  dp = dot((fn$normals), t(sn$normals))
  m$it[, dp < 0] = m$it[c(1, 3, 2), dp < 0]
  conv_sur <- m
  bbox <- meshcube(mesh)
  bbox_w <- sqrt(sum((bbox[1, ] - bbox[3, ])^2))
  bbox_h <- sqrt(sum((bbox[1, ] - bbox[2, ])^2))
  bbox_l <- sqrt(sum((bbox[1, ] - bbox[5, ])^2))
  bbox_c <- prod(bbox_w, bbox_h, bbox_l)
  voxelsize <- (bbox_c/ncells)^(1/3)
  xbox <- seq(min(bbox[, 1]), max(bbox[, 1]), by = voxelsize)
  ybox <- seq(min(bbox[, 2]), max(bbox[, 2]), by = voxelsize)
  zbox <- seq(min(bbox[, 3]), max(bbox[, 3]), by = voxelsize)
  GRID <- as.matrix(expand.grid(xbox, ybox, zbox))
  DISCR <- vcgClostKD(GRID, mesh, sign = TRUE)
  voxels_inside_1 <- GRID[which(DISCR$quality <= 0), ]
  out_endo <- vcgClostKD(voxels_inside_1, conv_sur, sign = TRUE)
  voxels_inside_2 <- voxels_inside_1[which(out_endo$quality >= 
                                             0), ]
  count_cells <- nrow(voxels_inside_2)
  vol <- ((voxelsize^3) * count_cells)/1000
  return(vol)
}

