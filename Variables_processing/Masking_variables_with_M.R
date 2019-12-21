#########
# Masking and saving reduced raster layers 
#########

# Description
## The following script helps to mask a set of raster layers to be used . 

## No data is needed if internet conection is available.

## Examples in which raster layers are masked to the polygons obtained (during) 
## the same process are also shown.

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


# assuming that you still don't have an M the occurrence data 
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