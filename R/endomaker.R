#' endomaker
#'
#' Build endocast from a skull 3D mesh 
#' @param mesh mesh3d: 3D model of the skull
#' @param path_in character: path of the skull where is stored
#' @param param1_endo numeric: parameter for spherical flipping 
#' @param npovs numeric: number of Points of View used in the endocast construction
#' @param decmesh numeric: decmesh
#' @param alpha_ext numeric: alpha shape for construction external cranial mesh
#' @param alpha_vol numeric: alpha shape for volume calculation 
#' @param plot logical: if TRUE the endocast is plotted
#' @param save logical: if TRUE the mesh of the endocast is saved
#' @param outpath character: path where save the endocast 
#' @param colmesh character: color of the mesh to be plotted
#' @param ncells numeric: approximative number of cell for 3D grid construction
#' @param npovs_calse numeric: number of Points of View for construction of skull shell
#' @param param1_calse numeric: parameter for calse (construction shell)
#' @param param1_ast numeric: parameter for ast3d (construction row endocast)
#' @param decendo numeric: desired number of triangles (row endocast)
#' @param scalendo numeric: scale factor row endocast (for definition of POVs)
#' @param alpha_end numeric: alpha shape value for concave hull (row endocast) 
#' @param mpovdist numeric: mean value between POVs and mesh 
#' @param volume logical: if TRUE the calculation of the volume (expressed in cc) through concave is returned
#' @param nVoxels numeric: number of voxels for estimation endocranial volume
#' @param num.cores numeric: numbers of cores to be used in parallel elaboration
#' @return endocast mesh3d: mesh of the endocast
#' @return volume numeric: volume of the endocast expressed in cc
#' @author Antonio Profico, Costantino Buzi, Marina Melchionna, Paolo Piras, Pasquale Raia, Alessio Veneziano
#' @references Profico, A., Buzi, C., Melchionna, M., Veneziano, A., & Raia, P. (2020). 
#' Endomaker, a new algorithm for fully automatic extraction of cranial endocasts and the calculation of their volumes. American Journal of Physical Anthropology.
#' @examples
#' \dontrun{
#' library(rgl)
#' data(human_skull)
#' sapendo<-endomaker(human_skull,param1_endo = 1.0,decmesh = 20000, num.cores=NULL)
#' open3d()
#' wire3d(sapendo$endocast,col="violet")
#' ecv<-sapendo$volume
#' }
#' @export

endomaker<-function (mesh = NULL, path_in = NULL, param1_endo = 1, npovs = 50, 
                     volume = TRUE, alpha_vol = 100, nVoxels = 1e+05, decmesh = 20000, 
                     alpha_ext = 30, ncells = 50000, npovs_calse = 50, param1_calse = 2, 
                     param1_ast = 1.3, decendo = 20000, scalendo = 0.5, alpha_end = 100, 
                     mpovdist = 10, plot = FALSE, colmesh = "orange", save = FALSE, 
                     outpath = tempdir(), num.cores = NULL) 
{


  if ((is.null(mesh) & is.null(path_in) == TRUE)) {
    stop
    cat("please, specify the mesh also input or path set where the mesh is stored")
  }
  if (is.null(mesh) == TRUE) {
    mesh <- vcgIsolated(file2mesh(path_in))
  }
  if (is.null(path_in) == TRUE) {
    mesh <- vcgIsolated(mesh)
  }
  dec_mesh <- vcgQEdecim(mesh, decmesh)
  ver_mesh <- vert2points(dec_mesh)
  conv_ver <- ashape3d(ver_mesh, alpha = alpha_ext, pert = TRUE, 
                       eps = 1e-09)
  x<-conv_ver
  selrows = x$triang[, 8 + 1] %in% 2L
  tr <- x$triang[selrows, c("tr1", "tr2", "tr3")]
  m=rgl::tmesh3d(
  vertices = t(x$x),
  indices = t(tr),
  homogeneous = FALSE)
  sn=alphashape3d::surfaceNormals(x, indexAlpha = 1)
  fn=facenormals(m)
  dp=dot((fn$normals), t(sn$normals))
  m$it[,dp<0]=m$it[c(1,3,2),dp<0]
  conv_sur<-m
  
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
  DISCR <- vcgClostKD(GRID, conv_sur, sign = TRUE)
  matrice <- GRID[which(DISCR$quality >= 0), ]
  DISCR2 <- vcgClostKD(GRID, mesh, sign = TRUE)
  matrice2 <- GRID[which(DISCR2$quality >= 0), ]
  field_mat <- vcgKDtree(matrice2, GRID, k = 5)
  dist_mat <- rowSums(field_mat$distance)
  QUANT <- quantile(dist_mat)
  pos_1 <- which(dist_mat < QUANT[2])
  pos_2 <- which(dist_mat >= QUANT[2] & dist_mat <= QUANT[3])
  pos_3 <- which(dist_mat > QUANT[3])
  cols <- NULL
  ca_lse <- ext.int.mesh(mesh = mesh, views = npovs_calse, 
                         param1 = param1_calse, default = TRUE, import_pov = NULL, 
                         num.cores = num.cores)
  vis <- out.inn.mesh(ca_lse, mesh, plot = FALSE)
  vis_mesh <- vis$visible
  GRID3 <- GRID[(pos_3), ]
  DISCR3 <- vcgClostKD(GRID3, vis_mesh, sign = TRUE)
  set_4 <- GRID3[which(DISCR3$quality >= 0), ]
  center_bone <- colMeans(GRID[pos_1, ])
  mat_center_bone <- repmat(center_bone, dim(set_4)[1], 3)
  mat_center_dists <- sqrt(rowSums((set_4 - mat_center_bone)^2))
  clusters <- kmeans(mat_center_dists, 2)
  a <- set_4[clusters$cluster == 1, ]
  b <- set_4[clusters$cluster == 2, ]
  ast3d <- ext.int.mesh(mesh = mesh, views = 0, param1 = param1_ast, 
                        default = FALSE, import_pov = TRUE, matrix_pov = t(as.matrix((colMeans(b)))), 
                        num.cores = num.cores)
  vis_inv_mesh <- out.inn.mesh(ast3d, mesh, plot = FALSE)
  vis_mesh <- (vis_inv_mesh$visible)
  raw_endo <- vcgQEdecim(vis_mesh, decendo)
  raw_endo_s <- scalemesh(raw_endo, scalendo, center = "m")
  raw_endo_v <- vert2points(raw_endo_s)
  conv_endo_1 <- ashape3d(raw_endo_v, alpha = alpha_end, pert = TRUE, 
                          eps = 1e-09)
  
  x<-conv_endo_1
  selrows = x$triang[, 8 + 1] %in% 2L
  tr <- x$triang[selrows, c("tr1", "tr2", "tr3")]
  m=rgl::tmesh3d(
    vertices = t(x$x),
    indices = t(tr),
    homogeneous = FALSE)
  sn=alphashape3d::surfaceNormals(x, indexAlpha = 1)
  fn=facenormals(m)
  dp=dot((fn$normals), t(sn$normals))
  m$it[,dp<0]=m$it[c(1,3,2),dp<0]
  conv_endo_s<-m
  
  rmesh = vcgClostKD(as.matrix(GRID), conv_endo_s, sign = TRUE)
  POVs <- as.matrix(kmeans(GRID[rmesh$quality >= 0, ], npovs, 
                           iter.max = 10000)$centers)
  size_ast <- mpovdist/mean(closemeshKD(POVs, mesh, k = 2)$quality)
  mesh_end_s <- scalemesh(mesh, size_ast, center = "m")
  povs_mean <- colMeans(vert2points(mesh))
  povs <- translate3d(POVs, -povs_mean[1], -povs_mean[2], -povs_mean[3])
  povs <- povs * size_ast
  povs <- translate3d(povs, povs_mean[1], povs_mean[2], povs_mean[3])
  ast3d2 <- ext.int.mesh(mesh = mesh_end_s, views = 0, param1 = param1_endo, 
                         default = FALSE, import_pov = TRUE, matrix_pov = povs, 
                         num.cores = num.cores)
  vis_inv_mesh2 <- out.inn.mesh(ast3d2, mesh_end_s, plot = FALSE)
  vis_mesh2 <- (vis_inv_mesh2$visible)
  meshmean <- colMeans(vert2points(mesh))
  endo_sur <- translate3d(vis_mesh2, -meshmean[1], -meshmean[2], 
                          -meshmean[3])
  endo_sur$vb[1:3, ] <- endo_sur$vb[1:3, ] * (1/size_ast)
  endo_sur <- translate3d(endo_sur, meshmean[1], meshmean[2], 
                          meshmean[3])
  endo_sur <- vcgIsolated(endo_sur)
  vol = NULL
  if (volume == TRUE) {
    vol <- volendo(endo_sur, alpha_vol = alpha_vol, ncells = ncells)
  }
  if (save == TRUE) {
    pathout <- outpath
    mesh2ply(endo_sur, "endocast")
  }
  if (plot == TRUE) {
    open3d()
    shade3d(endo_sur, col = colmesh)
  }
  closeAllConnections()
  out <- list(endocast = endo_sur, volume = vol)
}