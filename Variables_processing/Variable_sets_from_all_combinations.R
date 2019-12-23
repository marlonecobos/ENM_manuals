#####
# Creating sets of all potential combinations of a group of variables
#####

# Description
## The following script helps to create sets of environmental variables resulted
## from all potential combinations of all the ones originally considered. Sets 
## will have always more than two variables. This process has been recently 
## suggested as an option to test what set of variables among a group considered 
## to be appropriate for creating a model, results to be the best based on 
## further processes of model calibration.

## Note: All variables must have the same projection, extent, and resolution.
## We suggest to work with Geographic projections WGS84, with no planar projection.

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

# IF YOU HAVE THE DATA IN YOUR DIRECTORY AS DESCRIBED BELOw, USE THIS
# variables need to be saved in a subdirectory named "bio", variables must be in 
# ascii format (.asc)

## reading data
varaibles_list <- list.files(path = "bio", pattern = ".asc", # vector of variables
                             full.names = TRUE)

variables <- stack(varaibles_list) # stack of variables


# IF YOU DON'T HAVE THE DATA, USE THIS
## download data
variables <- getData("worldclim", var = "bio", res = 10)[[c(1, 10, 12, 16)]]

## crop the variables to an "invented" area for model calibration
ext <- extent(-105, -85, 5, 25)
variables <- crop(variables, ext)

## create the folder for saving cropped variables
dir.create("All_variables") 

## names of variables
variable_names <- paste0("All_variables/", names(variables), ".asc") 

## writing selected variables
for (i in 1:length(variable_names)) {
  writeRaster(variables[[i]], filename = variable_names[i], format = "ascii")
}


#######################################################################################
# Preparing sets of variables
##################

## functions help
help(kuenm_varcomb)

## preparing function's arguments
in_folder <- "All_variables" # name of folder with variables to be combined in distinct sets
out_folder <- "Var_combinations" # name of folder that will contain the sets 
                                 # of combined variables (the function creates this)
n_min <- 2 # minimum number of variables per combination
in_format <- "ascii" 
out_format <- "ascii"

## preparing all sets
combinations <- kuenm_varcomb(var.dir = in_folder, out.dir = out_folder, 
                              min.number = n_min, in.format = in_format,
                              out.format = out_format)