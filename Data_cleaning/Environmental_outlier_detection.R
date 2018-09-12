#####
# Detecting environmental outliers in species occurrences
#####

# data needed 
## (1) Species occurrences, in this example we will use a file named sp_occ.csv. 
## This file must contain the following columns (in that order): 
## Species_name, Longitud, Latitud.

## (2) Environmental raster layers, in this example we will use two ascii 
## layers located in the folder (bio) in our working directory. Our layers are 
## Mean annual temperature and Annual Precipitation These varaibles can be in 
## GeoTiff format as well.

## (3) The region of interest, accessioble area (M), or calibration area.
## This area can be a shapefile or a raster layer; for purposes of demonstration
## we will use a shapefile, however, lines of code for using a raster layer are
## written as comments.

## Note: Coordinates must coincide with the projection of the environmental data.
## Projection of accessible area and variables must also be the same.
## We suggest to work with Geographic projection WGS84, with no planar projection.


# loading needed packages (packages will be automatically installed if required)
pcakages <- c("raster", "rgdal", "rgeos")
req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
if (length(req_packages) > 0) {
  install.packages(req_packages, dependencies = TRUE)
}
sapply(pcakages, require, character.only = TRUE)

# defining working directory
setwd("D:/Marlon/Data_cleaning/") # change this to your working directory

# reading data
occurrences <- read.csv("sp_occ.csv") # occurrences

variables <- stack("bio/bio_1", "bio/bio_12") # stack of variables

M <- readOGR(dsn = ".", layer = "acc_area")# accessible area (M), or calibration area

# masking variables to the calibration area
m_variables <- mask(varaibles, M)

# getting data from the variables
variables_values <- na.omit(values(m_variables)) # for the region of interest

occ_variables <- na.omit(extract(m_variables, occ[, 2:3]))

# sample of 10000 values if more pixels exist
if (dim(variables_values > 10000)) {
  variables_values <- variables_values[sample(rownames(variables_values), 10000), ] 
} 

# plot for searching for potential environmental outliers
par(mar = c(4.5, 4, 0.5, 0.5), cex = 0.9)
plot(variables_values[, 1], variables_values[, 2], col = "gray75",
     pch = 1, xlab = "Temperature", ylab = "Precipitation")
points(occ_variables[, 1], occ_variables[, 2], col = "black",
       pch = 19)
legend("topleft", legend = c("Region of interest", "Occurrences"),
       pch = c(1, 19), col = c("gray75", "black"))

# if need to save the figure use the following code
jpeg(filename = "Environmental_outliers.jpg", width = 166, height = 166,
     units = "mm", res = 600)
par(mar = c(4.5, 4, 0.5, 0.5), cex = 0.9)
plot(variables_values[, 1], variables_values[, 2], col = "gray75",
     pch = 1, xlab = "Temperature", ylab = "Precipitation")
points(occ_variables[, 1], occ_variables[, 2], col = "black",
       pch = 19)
legend("topleft", legend = c("Region of interest", "Occurrences"),
       pch = c(1, 19), col = c("gray75", "black"))
dev.off()  