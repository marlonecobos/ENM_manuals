#####
# Detecting environmental outliers in species occurrences
#####

# data needed 
## (1) Species occurrences; in this example we will use a file named sp_occ.csv. 
## This file must contain the following columns (in that order): 
## ID, Species_name, Longitud, Latitud.

## (2) Environmental raster layers; in this example we will use two ascii 
## layers located in the folder (bio) in our working directory. Our layers are 
## Mean annual temperature and Annual Precipitation. These varaibles can be in 
## GeoTiff format as well.

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
setwd("D:/Marlon/R/Data_cleaning") # change this to your working directory

# reading data
occurrences <- read.csv("sp_occ.csv") # occurrences

variables <- stack("bio/bio_1.asc", "bio/bio_12.asc") # stack of variables

M <- readOGR(dsn = ".", layer = "acc_area")# accessible area (M), or calibration area

# if M is a raster layer use
# M <- raster("acc_area.asc")

# masking variables to the calibration area
m_variables <- crop(variables, M)
m_variables <- mask(m_variables, M)

# getting data from the variables
variables_values <- na.omit(values(m_variables)) # for the region of interest

occ_variables <- na.omit(cbind(occurrences[, 1], # for occurrences adding IDs
                               extract(m_variables, occurrences[, 3:4])))

# sample of 10000 values if more pixels exist
if (dim(variables_values)[1] > 10000) {
  variables_values <- variables_values[sample(1:dim(variables_values)[1], 10000), ] 
} 

# plot for searching for potential environmental outliers
par(mfcol = c(1, 2), mar = c(4.5, 4, 0.5, 0.5), cex = 0.9)
plot(variables_values[, 1], variables_values[, 2], col = "gray65",
     pch = 1, xlab = "Temperature", ylab = "Precipitation")
points(occ_variables[, 2], occ_variables[, 3], col = "black",
       pch = 19)
legend("topleft", legend = c("Region of interest", "Occurrences"),
       pch = c(1, 19), col = c("gray65", "black"), bty = "n")
plot(variables_values[, 1], variables_values[, 2], col = "gray65",
     pch = 1, xlab = "Temperature", ylab = "")
text(occ_variables[, 2], occ_variables[, 3], occ_variables[, 1], cex = 0.6)
legend("topleft", legend = "Occurrence ID", bty = "n")

# if need to save the figure use the following code
jpeg(filename = "Environmental_outliers.jpg", width = 166, height = 83,
     units = "mm", res = 600)
par(mfcol = c(1, 2), mar = c(4.5, 4, 0.5, 0.5), cex = 0.65)
plot(variables_values[, 1], variables_values[, 2], col = "gray65",
     pch = 1, xlab = "Temperature", ylab = "Precipitation")
points(occ_variables[, 2], occ_variables[, 3], col = "black",
       pch = 19)
legend("topleft", legend = c("Region of interest", "Occurrences"),
       pch = c(1, 19), col = c("gray65", "black"), bty = "n")
plot(variables_values[, 1], variables_values[, 2], col = "gray65",
     pch = 1, xlab = "Temperature", ylab = "")
text(occ_variables[, 2], occ_variables[, 3], occ_variables[, 1], cex = 0.6)
legend("topleft", legend = "Occurrence ID", bty = "n")
dev.off()  