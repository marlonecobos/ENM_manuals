#' Principal componens for raster layers and projections
#' 
#' @description kuenm_rpca performs a principal component analysis with a set of variables and
#' produces raster layers of them. If needed the pricipal components are projected to other 
#' scenarios.
#' 
#' @param vars.folder (character) name of the folder where the environmental variables are. Variables
#' must be in ascii format and at least two.
#' @param in.format (character) format of variables in \code{var.dir}. Options are "ascii", "GTiff", and "EHdr" = bil.
#' Default = "ascii".
#' @param out.format (character) format of variables to be written in distinct sets inside \code{out.dir}.
#' Options are "ascii", "GTiff", and "EHdr" = bil. Default = "ascii". 
#' @param project (logical) whether or not to project the species niche to other scenario(s).
#' If TRUE, argument \code{proj.variables} needs to be defined. Default = FALSE.
#' @param multiple.proj (logical) whether or not the projection should be performed for various
#' scenarios. If TRUE, argument \code{proj.vars} changes; see below. Default = FALSE.
#' @param proj.vars (character) If \code{multiple.proj} is FALSE, name of the folder where  
#' environmental variables of the scenario for projection are. If \code{multiple.proj} is TRUE, name of 
#' the folder where folders with environmental variables of the scenarios for projections are. Variables 
#' must be in ascii format and they must correspond with variables in \code{vars.folder} (i.e., their name 
#' must correspond but they should represent conditions in the LGM).
#' @param return.in (logical) whether or not return raster layers of principal components to
#' the R environment in a list with other results.
#' @param n.pcs (numeric) if \code{proj.variables} = TRUE, number of principal components to be 
#' returned as rasters. By default all principal components are returned as layers.
#' @param out.dir (character) name of the folder to be created to save the results of the analyses.
#' Default = "PCA_results".
#' 
#' @return 
#' A list containing the PCA loading and PCA summary as a matrix; if \code{return.in} = TRUE a RasterLayer
#' (if not projected) or a RasterStack (if projected) of the principal components are returned additionally.
#' All results are written in \code{out.dir}.
#' 
#' @details 
#' 
#' 
#' 
#' @examples 
#' # Arguments


kuenm_rpca <- function(vars.folder, in.format = "ascii", out.format = "ascii", project = FALSE, 
                       multiple.proj = FALSE, proj.vars, return.in = FALSE, n.pcs, out.dir = "PCA_results") {
  pcakages <- c("raster")
  req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
  if (length(req_packages) > 0) {
    install.packages(req_packages, dependencies = TRUE)
  }
  
  if (missing(vars.folder)) {
    stop("Argument vars.folder must be defined. See functions help.")
  }
  if (project == TRUE) {
    if (missing(proj.vars)) {
      stop("If projections are needed, argument proj.vars must be defined. See functions help.")
    }
  }
  if (in.format == "ascii") {
    patt <- ".asc$"
  }
  if (in.format == "GTiff") {
    patt <- ".tif$"
  }
  if (in.format == "EHdr") {
    patt <- ".bil$"
  }
  if (out.format == "ascii") {
    patt1 <- ".asc"
  }
  if (out.format == "GTiff") {
    patt1 <- ".tif"
  }
  if (out.format == "EHdr") {
    patt1 <- ".bil"
  }
  
  # reading variables
  var <- list.files(vars.folder, pattern = patt, full.names = TRUE)
  
  variab <- raster::stack(var)
  var_points <- na.omit(raster::values(variab))
  
  # pca analyses
  pca <- prcomp(var_points, center = TRUE, scale = TRUE)
  scores <- pca$x
  
  dir.create(out.dir)
  
  pca_fol <- paste(out.dir, "Initial", sep = "/")
  dir.create(pca_fol)
  
  if (return.in == TRUE) {
    pcras <- list()
    if (missing(n.pcs)) {
      n.pcs <- length(var)
    }
  }
  
  for (i in 1:n.pcs) {
    pcra <- variab[[1]]
    pcra[!is.na(raster::values(pcra))] <- scores[, i]
    filenam <- paste(pca_fol, "/pc_", i, patt1, sep = "")
    raster::writeRaster(pcra, filenam, format = out.format)
    
    if (return.in == TRUE) {
      pcras[[i]] <- pcra
    }
  }
  
  if (return.in == TRUE) {
    pcras <- do.call(raster::stack, pcras)
  }
  
  StdDev <- pca$sdev
  VarExp <- pca$sdev^2/sum(pca$sdev^2)
  CumVar <- cumsum(VarExp)
  SumPCAMat <- rbind(StdDev, VarExp, CumVar)
  colnames(SumPCAMat) <- paste("PC", seq(1, length(StdDev)), sep = "")
  row.names(SumPCAMat) <- c("Standard deviation", "Proportion of Variance", 
                            "Cumulative Proportion")
  
  sink(paste(paste(pca_fol, "pca_results.txt", sep = "/")))
  cat("Principal component analysis results\n")
  cat("\nPCA loadings\n")
  print(pca$rotation)
  
  cat("\n\nPCA summary\n")
  print(SumPCAMat)
  sink()
  
  # pca results to be returned
  loadings <- pca$rotation
  respca <- SumPCAMat
  
  # projecting PCs
  if (project == TRUE) {
    proj_dirs <- list.dirs(proj.vars, recursive = FALSE)
    proj_names <- list.dirs(proj.vars, recursive = FALSE, full.names = FALSE)
    fol_names <- paste(out.dir, proj_names, sep = "/")
    
    if (return.in == TRUE) {
      ppcrass <- list()
    }
    
    for (h in 1:length(proj_dirs)) {
      pvar <- list.files(proj_dirs[h], pattern = patt, full.names = TRUE)
      p_stack <- raster::stack(pvar)
      dir.create(fol_names[h])
      
      if (return.in == TRUE) {
        ppcras <- list()
      }
      
      p_stackp <- na.omit(raster::values(p_stack))
      names(p_stackp) <- names(pca[[4]])
      p_pcs <- predict(pca, newdata = p_stackp)
      
      for (i in 1:n.pcs) {
        pcra <- p_stack[[1]]
        pcra[!is.na(raster::values(pcra))] <- p_pcs[, i]
        filenam <- paste(fol_names, "/pc_", i, patt1, sep = "")
        raster::writeRaster(pcra, filenam, format = out.format)
        
        if (return.in == TRUE) {
          ppcras[[i]] <- pcra
        }
      }
      
      if (return.in == TRUE) {
        ppcrass[[h]] <- do.call(raster::stack, ppcras)
      }
    }
    names(ppcrass) <- paste("PCRasters", proj_names, sep = "_")
  }
  
  if (return.in == TRUE) {
    if (project == TRUE) {
      results <- c(list(loadings, respca, pcras), ppcrass)
      names(results)[1:3] <- c("PCA_loadings", "PCA_results", "PCRasters_initial")
    }else {
      results <- list(loadings, respca, pcras)
      names(results) <- c("PCA_loadings", "PCA_results", "PCRasters_initial")
    }
  }else {
    results <- list(loadings, respca)
    names(results) <- c("PCA_loadings", "PCA_results")
  }
  
  return(results)
}
