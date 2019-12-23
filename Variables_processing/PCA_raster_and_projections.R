#####
# Creating raster pcs and their projections (if needed) from raster variables
#####

# Description
## The following script helps to perform a Pricipal Component Analysis (PCA) on
## raster layers. The product are statistical results and raster Principal 
## Components that can be used as predictors in further processes. This script 
## also helps in projecting such Principal Components to other raster layers,
## so the Principal Components for the projected layers can be comparable to the 
## initial ones.

## Note: Variables from each scenario must have the same projection, extent, and 
## resolution. We suggest to work with Geographic projections WGS84, with no 
## planar projection.

## The main processes are performed the package kuenm from GitHub. To install
## this package see instructions in https://github.com/marlonecobos/kuenm.

# loading needed packages (packages will be automatically installed if required)
suppressWarnings({
  if(!require(raster)){
    install.packages("raster")
    library(raster)
  }
}) 

# loading needed package
# assuming that you installed kuenm, load it, if not installed see 
# https://github.com/marlonecobos/kuenm for instructions
library(kuenm)


#######################################################################################
# Preparing directory and data
##################

# defining working directory
setwd("D:/Marlon/Variables_processing/") # change this to your working directory

# IF YOU HAVE THE DATA IN YOUR DIRECTORY AS DESCRIBED BELOW, USE THIS
# variables need to be saved in a subdirectory named "bio", variables must be in 
# ascii format (.asc)

## reading data
varaibles_list <- list.files(path = "bio", pattern = ".asc", # vector of variables
                             full.names = TRUE)

variables <- stack(varaibles_list) # stack of variables


# IF YOU DON'T HAVE THE DATA, USE THIS
## download data
variables <- getData("worldclim", var = "bio", res = 10)[[-c(8, 9, 18, 19)]]

## crop the variables to an "invented" area for model calibration
ext <- extent(-105, -85, 5, 25)
cvariables <- crop(variables, ext)

## create the folder for saving cropped variables
dir.create("Crop_variables") 

## names of variables
variable_names <- paste0("Crop_variables/", names(cvariables), ".asc") 

## writing selected variables
for (i in 1:length(variable_names)) {
  writeRaster(cvariables[[i]], filename = variable_names[i], format = "ascii")
}

# for projections
## create the folder for saving world variables
dir.create("World_variables") 

## names of variables
variable_names <- paste0("World_variables/", names(variables), ".asc") 

## writing selected variables
for (i in 1:length(variable_names)) {
  writeRaster(variables[[i]], filename = variable_names[i], format = "ascii")
}


#######################################################################################
# Preparing sets of variables
##################

# simple raster PCA
## functions help
help(kuenm_rpca)

## preparing function's arguments
var_folder <- "Crop_variables" # name of folder with variables to be combined in distinct sets
in_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil 
scalev <- TRUE # scale variables
writer <- TRUE # save results
out_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
out_folder <- "PCA_results" # name of folder that will contain the sets 
n_pcs <- 6 # number of pcs you want as rasters, if not defined all pcs are returned as rasters


## runing PCA
kuenm_rpca(variables = var_folder, in.format = in_format, var.scale = scalev, 
           write.result = writer, out.format = out_format, out.dir = out_folder,
           n.pcs = n_pcs)



# raster PCA with projections
## preparing function's arguments
var_folder <- "Crop_variables" # name of folder with variables to be combined in distinct sets
in_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil 
scalev <- TRUE # scale variables
writer <- TRUE # save results
out_format <- "ascii" # other options available are "GTiff" and "EHdr" = bil
out_folder <- "PCA_results_proj" # name of folder that will contain the sets 
project <- TRUE
proj_folder <- "World_variables"
n_pcs <- 6 # number of pcs you want as rasters, if not defined all pcs are returned as rasters


## runing PCA
kuenm_rpca(variables = var_folder, in.format = in_format, var.scale = scalev, 
           write.result = writer, out.format = out_format, out.dir = out_folder,
           project = project, proj.vars = proj_folder, n.pcs = n_pcs)