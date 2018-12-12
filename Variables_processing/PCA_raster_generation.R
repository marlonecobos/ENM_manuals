#####
# Creating raster pcs and their projections (if needed) from raster variables
#####

# data needed 
## Environmental raster layers, in this example we will use a set of ascii 
## layers located in the folder (reduced_variables) in our working directory. 

## If pcs need to be projected to other scenarios, then a folder containing folders with 
## variables representing other scenarios is needed as well.

## Variables can be ascii, tif, or bil files.

## Note: All variables must have the same projection, extent, and resolution.
## We suggest to work with Geographic projections WGS84, with no planar projection.

source("https://raw.githubusercontent.com/marlonecobos/ENM_manuals/master/Variables_processing/kuenm_rpca.R")

setwd("YOUR/DIRECTORY")

# no variable projections to distinct scenarios needed

var_folder <- "reduced_variables" # name of folder with variables to be combined in distinct sets
out_folder <- "PCA_results" # name of folder that will contain the sets 
in_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil 
out_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
n_pcs <- 6 # number of pcs you want as rasters, if not defined all pcs are returned as rasters

kuenm_rpca(vars.folder = var_folder, in.format = "ascii", out.format = "ascii", project = FALSE, 
           n.pcs = n_pcs, out.dir = out_folder)



# ifvariable projections to distinct scenarios are needed

var_folder <- "reduced_variables" # name of folder with variables to be combined in distinct sets
proj_folder <- "scenarios" # name of the folder containing one or more folders with variables for other scenarios
out_folder <- "PCA_results" # name of folder that will contain the sets 
in_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil 
out_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
n_pcs <- 6 # number of pcs you want as rasters, if not defined all pcs are returned as rasters

kuenm_rpca(vars.folder = var_folder, in.format = "ascii", out.format = "ascii", project = TRUE, 
           proj.vars = proj_folder, n.pcs = n_pcs, out.dir = out_folder)