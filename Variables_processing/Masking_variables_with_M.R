#########
# Masking and saving reduced raster layers 
#########

# Description
## The following script helps to mask a set of raster layers to be used in further
## analyses in Ecological niche modeling exercises. 

## No data is needed if internet conection is available.

## The main processes are performed the package ellisenm from GitHub. To install
## this package see instructions in https://github.com/marlonecobos/ellipsenm.

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
  if(!require(raster)){
    install.packages("raster")
    library(rgbif)
  }
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

# raster layes (only 6 coarse resolution)
variables <- getData("worldclim", var = "bio", res = "10")[[c(1, 10, 11, 12, 16, 17)]] # stack of variables


# assuming that you still don't have an M, here is an example 
sp <- "Cynomys mexicanus" # species name

occ <- occ(query = sp, from = "gbif", limit = 1000) # getting data
occ <- fixnames(occ, how = "query")$gbif$data[[1]] # fix_names

## keeping only unique georeferenced records.
occ_g <- occ[!is.na(occ$longitude) & !is.na(occ$latitude), ] # excluding no georeferences
occ_g <- occ_g[!duplicated(paste(occ_g$name, occ_g$longitude, # excluding duplicates
                                 occ_g$latitude, sep = "_")), ]
occ_g <- occg[, c("name", "longitude", "latitude")] # only these three columns

## the following line will spatially rarefy the data (thinning) distance = 10 km
occ_t <- thin_data(occ_g, longitude = "longitude", latitude = "latitude", 
                   thin_distance = 10)

## the M
M <- convex_area(occ_t, longitude = "longitude", latitude = "latitude", 
                 buffer_distance = 75)


#######################################################################################
# Masking and saving layers
##################

# masking variables
var_mask <- mask(crop(variables, M_intersect), M_intersect)

# saving masked variables  
## new directory
dir.create("Masked_layers")

## names for layers
rnames <- paste0("Masked_layers/", names(variables), ".asc") # users select the format

## saving layers in new folder
sav <- lapply(1:nlayers(var_mask), function(x) {
  writeRaster(var_mask[[x]], filename = rnames, format = "ascii") # change format accordingly
})