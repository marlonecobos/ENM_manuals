#########
# Creating hypothesis of areas for model calibration (M)
#########

# Description
## The following script helps to create spatial polygons to represent areas for
## calibration of ecological niche models four types of areas will be created in
## this example: (1) using buffers; (2) creating convex hulls, (3) creating
## concave hulls, and (4) selecting polygons from a layer of spatial plygons that
## represent ecorregions. 

## No data is needed if internet conection is available.

## One of the main processes is performed the package ellisenm from GitHub. To 
## install it see instructions in https://github.com/marlonecobos/ellipsenm.


# loading needed package
suppressWarnings({
  if(!require(spocc)){
    install.packages("spocc")
    library(spocc)
  }
  #if(!require(scrubr)){
  #  install.packages("scrubr")
  #  library(scrubr)
  #}
}) 

# assuming that you installed ellipsenm, load it, if not installed see 
# https://github.com/marlonecobos/ellipsenm for instructions
library(ellipsenm)


#######################################################################################
# Preparing directory and data
##################

# defining working directory
## project forlder
setwd("Z:/Marlon_E_Cobos/ENM_project") # Your folder


# getting the occurrence data 
sp <- "Cynomys mexicanus" # species name

occ <- occ(query = sp, from = "gbif", limit = 1000) # getting data
occ <- fixnames(occ, how = "query")$gbif$data[[1]] # fix_names

## keeping only unique georeferenced records.
occ_g <- occ[!is.na(occ$longitude) & !is.na(occ$latitude), ] # excluding no georeferences
occ_g <- occ_g[!duplicated(paste(occ_g$name, occ_g$longitude, # excluding duplicates
                                 occ_g$latitude, sep = "_")), ]
occ_g <- occg[, c("name", "longitude", "latitude")] # only these three columns


#######################################################################################
# Preparing directory and data
##################

# The results of these processes can be saved to the working directory directly.
# To do so, check the arguments "save" and "name". 

# Thinning (spatially rarefaction of data distance = 10 km)
occ_t <- thin_data(occ_g, longitude = "longitude", latitude = "latitude", 
                   thin_distance = 10)


# Splitting data into training and testing
## This function will help splitting records to perform processes of model 
## calibration, a list with all, train, and test records will be returned
occ_calibration <- ellipsenm::split_data(occ_t, method = "random", 
                                         longitude = "longitude", latitude = "latitude", 
                                         train_proportion = 0.75)