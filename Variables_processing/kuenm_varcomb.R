#' All potential combinations of a group of variables
#'
#' @description kuenm_varcomb creates multiple sets of variables by grouping them in all their potential combinations.
#'
#' @param var.dir (character) the name of the folder containing variables that will be combined.
#' @param out.dir (character) the name of the folder in which subfolders with distinct combinations of
#' variables will be written.
#' @param min.number (integer) the minimum number of variables per combination. This number must be > 1. 
#' Default = 2.
#' @param in.format (character) format of variables in \code{var.dir}. Options are "ascii" and "GTiff". 
#' Default = "ascii".
#' @param out.format (character) format of variables to be written in distinct sets inside \code{out.dir}. 
#' Options are "ascii" and "GTiff". Default = "ascii".
#'
#' @return A list containing vectors of all the potential combinations of variables. In addition, a folder 
#' named \code{out.dir} with subfolders in which distinct combinations of variables produced are written.
#'
#' @details Sest of variables are written in the working directory and not retained as RasterStacks to avoid
#' problems related to RAM limitations.


kuenm_varcomb <- function(var.dir, out.dir, min.number = 2, in.format = "ascii", 
                          out.format = "ascii") {
  
  # Setting things up
  if (min.number < 2) {
    stop("min.number must be an integer > 1.")
  }
  if (in.format == "ascii") {
    patt <- ".asc$"
  }
  if (in.format == "GTiff") {
    patt <- ".tif$"
  }
  if (out.format == "ascii") {
    patt1 <- ".asc$"
  }
  if (out.format == "GTiff") {
    patt1 <- ".tif$"
  }
  
  # List variable names
  variables <- list.files(path = var.dir, pattern = patt)
  
  if (length(variables) == 0) {
    stop(paste("No variables with format", in.format, "were found in the directory", var.dir))
  }
  
  if (min.number > length(variables)) {
    stop("min.number must be < the total number of variables.")
  }
  
  # Generating all combinations of variable names
  var_comb <- list()
  
  for (i in min.number:length(variables)) {
    comb <- combn(x = variables, m = i)
    comb_vs <- list()
    
    for (j in 1:dim(comb)[2]) {
      comb_vs[[j]] <- comb[, j]
    }
    
    var_comb[[i]] <- comb_vs
  }
  
  var_comb[1:(min.number - 1)] <- NULL
  var_combinations <- do.call(c, var_comb)
  
  # Preparing folders, variable combinations, and writing results
  dir.create(out.dir)
  
  vars_all <- raster::stack(paste(var.dir, variables, sep = "/"))
  
  sub_paths <- paste(out.dir, paste("Set", 1:length(var_combinations), sep = "_"), sep = "/")
  
  cat("\nA total of", length(sub_paths), "sets of variables resulted from combinations of", length(variables), "variables will be written.\n")
  
  if(.Platform$OS.type == "unix") {
    pb <- txtProgressBar(min = 0, max = length(sub_paths), style = 3) #progress bar
  } else {
    pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sub_paths), width = 300) #progress bar
  }
  
  for (k in 1:length(sub_paths)) {
    Sys.sleep(0.1)
    if(.Platform$OS.type == "unix") {
      setTxtProgressBar(pb, i)
    } else {
      setWinProgressBar(pb, k, title = paste(round(k / length(sub_paths) * 100, 2),
                                             paste("% of the process has finished")))
    }
    
    dir.create(sub_paths[k])
    vars_set <- vars_all[[gsub(patt, "", var_combinations[[k]])]]
    
    for (l in 1:dim(vars_set)[3]) {
      raster::writeRaster(vars_set[[l]], filename = paste(sub_paths[k], var_combinations[[k]][l], sep = "/"), 
                          format = out.format)
    }
  }
  
  if(.Platform$OS.type != "unix") {
    suppressMessages(close(pb))
  }
  
  return(var_combinations)
}


