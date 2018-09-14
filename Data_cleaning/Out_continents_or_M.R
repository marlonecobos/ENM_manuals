#####
# Detecting occurrences outside the continents or the region of interest
#####

# data needed 
## (1) Species occurrences; in this example we will use a file named sp_occ.csv. 
## This file must contain the following columns (in that order): 
## ID, Species_name, Longitud, Latitud.

## (2) Layer representing the continents; in this example we will use an ascii 
## layer located in the folder (bio) in our working directory. Our layer is 
## Mean annual temperature. This varaible can be in GeoTiff format as well.
## We recommend to use a raster layer better, as this layer could be one of 
## the set of varaibles that will be used when creating the models.

## (3) The region of interest, accessible area (M), or calibration area.
## This area can be a shapefile or a raster layer; for purposes of demonstration
## we will use a shapefile, however, lines of code for using a raster layer are
## written as comments.

## Note: Coordinates must coincide with the projection of the environmental data.
## Projection of accessible area and variables must also be the same projection.
## We suggest to work with Geographic projection WGS84, with no planar projection.


# loading needed packages (packages will be automatically installed if required)
pcakages <- c("raster", "rgdal")
req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
if (length(req_packages) > 0) {
  install.packages(req_packages, dependencies = TRUE)
}
sapply(pcakages, require, character.only = TRUE)

# defining working directory
setwd("C:/Users/Marlon/Documents/R/Data_cleaning") # change this to your working directory

# reading data
occurrences <- read.csv("Cynomys_ludovicianus_Exercise.csv") # occurrences

continents <- raster("bio/bio_1.asc") # representing the continents

m <- readOGR(dsn = ".", layer = "acc_area") # accessible area (M), or calibration area
M <- crop(rasterize(m, continents), m)

# if M is a raster layer use
# M <- raster("acc_area.asc")

# checking which records are out of the continents
occ_inout <- data.frame(occurrences, 
                        inside = extract(continents, occurrences[, 3:4]))
occ_inout$inside <- ifelse(!is.na(occ_inout$inside), # only TRUE are inside
                           TRUE, FALSE)

## plotting will help to decide which ones to move
par(mar = c(0.5, 0.5, 0.5, 0.5), cex = 0.9)
image(continents)
box()
points(occ_inout[, 3], occ_inout[, 4], col = "black",
       pch = 19)

## excluding all records outside the continent
occurrences <- occurrences[occ_inout$inside == TRUE, ]

# checking which records are out of the region of interest
occ_reg_inout <- data.frame(occurrences, 
                        inside = extract(M, occurrences[, 3:4]))
occ_reg_inout$inside <- ifelse(!is.na(occ_reg_inout$inside), # only TRUE are inside
                           TRUE, FALSE)

## plotting will help to decide which ones to move
par(mar = c(0.5, 0.5, 0.5, 0.5), cex = 0.9)
image(continents)
image(M, col = "blue", add = TRUE)
box()
points(occ_reg_inout[, 3], occ_reg_inout[, 4], col = "black",
       pch = 19)

## excluding all records outside the continent
occurrences <- occurrences[occ_reg_inout$inside == TRUE, ]

# saving the new set of occurrences inide continents and area of interest
write.csv(occurrences, "Cynomys_ludovicianus_Exercise_inside.csv", row.names = FALSE)