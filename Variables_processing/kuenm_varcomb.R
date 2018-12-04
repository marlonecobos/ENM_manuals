#' All potential combinations of a group of variables
#'
#' @description kuenm_varcomb creates multiple sets of variables by grouping them in all their potential combinations.
#'
#' @param vars (character) the name of the folder containing variables (.asc format) that will be combined.
#' @param out.vars (character) the name of the folder in which subfolders with distinct combinations of
#' variables (.asc format) will be written.
#' @param min.number (integer) the minimum number of variables per combination. This number must be > 1. 
#' Default = 2.
#'
#' @return A list containing vectors of all the potential combinations of vars. In addition, a folder named out.vars 
#' with subfolders in which the distinct combinations of variables produced were written.
#'
#' @details Sest of variables are written in the working directory and not retained as RasterStacks to avoid
#' problems related to RAM limitations.


kuenm_varcomb <- function(vars, out.vars, min.number = 2) {
  
  if (min.number < 2) {
    stop("min.number must be an integer > 1.")
  }
  
  # List variable names
  variables <- list.files(path = paste(getwd(), vars, sep = "/"), pattern = ".asc")
  
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
  path <- paste(getwd(), out.vars, sep = "/")
  dir.create(path)
  
  vars_all <- raster::stack(paste(getwd(), vars, variables, sep = "/"))
  
  sub_paths <- paste(path, paste("Set", 1:length(var_combinations), sep = "_"), sep = "/")
  
  pb <- winProgressBar(title = "Progress bar", min = 0, max = length(sub_paths), width = 300) #progress bar
  
  for (k in 1:length(sub_paths)) {
    Sys.sleep(0.1)
    setWinProgressBar(pb, k, title = paste(round(k / length(sub_paths) * 100, 2),
                                           paste("% of the process has finished")))
    
    dir.create(sub_paths[k])
    vars_set <- vars_all[[gsub(".asc", "", var_combinations[[k]])]]
    
    for (l in 1:dim(vars_set)[3]) {
      raster::writeRaster(vars_set[[l]], filename = paste(sub_paths[k], var_combinations[[k]][l], sep = "/"), 
                          format = "ascii")
    }
  }
  suppressMessages(close(pb))
  
  return(var_combinations)
}


