#########
#  M by intersection of multiple M-hypotheses
#########

# Description
## The following script helps to create an area for model calibration based on 
## the intersection of multiple hypotheses of M. 

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
  if(!require(rgdal)){
    install.packages("rgdal")
    library(rgdal)
  }
  if(!require(rgeos)){
    install.packages("rgeos")
    library(rgeos)
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


# getting the occurrence data 
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

# getting ecorregions
download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
              destfile = file.path(getwd(), "wwf_ecoregions.zip")) 
unzip(file.path(getwd(), "wwf_ecoregions.zip"), 
      exdir = file.path(getwd(), "WWF_ecoregions"))
file.remove(file.path(getwd(), "wwf_ecoregions.zip"))
ecor <- readOGR("WWF_ecoregions/official", layer = "wwf_terr_ecos")

#######################################################################################
# Creating Ms as initial hypotheses
##################

# Before using the following functions, make sure you check their documentation
# (healp). They have handly options that allow you to mask raster layers as well 
# as to save the results directly. 

# Hint: use argument "save" to write the results in your directory, if needed

#####
# Areas by buffering records (100 km buffer)
M_buffer <- buffer_area(occ_t, longitude = "longitude", latitude = "latitude", 
                        buffer_distance = 100)

#####
# Areas using convex hulls (including 75 km buffer)
M_convex <- convex_area(occ_t, longitude = "longitude", latitude = "latitude", 
                        buffer_distance = 75)

#####
# Areas by selecting polygons (including 25 km buffer)
M_ecorreg <- polygon_selection(occ_t, longitude = "longitude", latitude = "latitude",
                               polygons = ecor, buffer_distance = 25)


#######################################################################################
# Intersecting to obtain an consensus M
##################

# intersection
M_intersect <- gIntersection(M_buffer, M_convex)
M_intersect <- gIntersection(M_intersect, M_ecorreg)

# visualization
par(mfrow = c(2, 2), cex = 0.6, mar = rep(0.3, 4))
plot(M_buffer); points(occ_t[, 2:3]); legend("topleft", legend = "Buffer", bty = "n")
plot(M_convex); points(occ_t[, 2:3]); legend("topleft", legend = "Convex hull", bty = "n")
plot(M_ecorreg); points(occ_t[, 2:3]); legend("topleft", legend = "Ecorregions", bty = "n")
plot(M_intersect); points(occ_t[, 2:3]); legend("topleft", legend = "Intersection", bty = "n")


# saving final product
writeOGR(M_intersect, ".", "M_intersection", driver = "ESRI Shapefile")
