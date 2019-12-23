#####
# Evaluation of variables correlation
#####

# Description
## The following script helps to measure correlation among distinct raster layers
## to be used as predictors in ecological niche modeling. This process needs to
## be performed in the area of model calibration.

## Note: All variables must have the same projection, extent, and resolution.
## We suggest to work with Geographic projections WGS84, with no planar projection.

## The main processes are performed the package ellisenm from GitHub. To install
## this package see instructions in https://github.com/marlonecobos/ellipsenm.

# loading needed packages (packages will be automatically installed if required)
suppressWarnings({
  if(!require(raster)){
    install.packages("raster")
    library(raster)
  }
}) 

# assuming that you installed ellipsenm, load it, if not installed see 
# https://github.com/marlonecobos/ellipsenm for instructions
library(ellipsenm)


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
variables <- crop(variables, ext)


#######################################################################################
# Varriable correlation analysis
##################

# functions help
help(variable_correlation)
# to save matrix correlation results see the arguments "save" and "name"

# correlation matrix
cors <- variable_correlation(variables)

# checking the table
View(cors) 

# analysis and ploting (all values above 0.75 will be magnified)
cors1 <- variable_correlation(variables, correlation_limit = 0.75,
                              corrplot = TRUE, magnify_to = 1.5)


#######################################################################################
# Selecting variables
##################

# selecting variables and writing them in a new directory
names(variables) # names

selected <- c("bio2", "bio4", "bio7", "bio13", "bio17") # change to the ones you select

selected_variables <- variables[[selected]] # selecting

## create the folder for selected variables
dir.create("Non_correlated_variables") 

## names of variables
variable_names <- paste0("Non_correlated_variables/", selected, ".asc") 

## writing selected variables
for (i in 1:nlayers(selected_variables)) {
  writeRaster(selected_variables[[i]], filename = variable_names[i], format = "ascii")
}
