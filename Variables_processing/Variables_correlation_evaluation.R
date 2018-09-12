#####
# Detecting environmental outliers in species occurrences
#####

# data needed 
## Environmental raster layers, in this example we will use a set of ascii 
## layers located in the folder (bio) in our working directory. Our layers are 
## 10 bioclimatic varaibles. These varaibles can be in GeoTiff format as well.

## Note: All variables must have the same projection.
## We suggest to work with Geographic projections WGS84, with no planar projection.


# loading needed packages (packages will be automatically installed if required)
pcakages <- c("raster", "rgdal")
req_packages <- pcakages[!(pcakages %in% installed.packages()[, "Package"])]
if (length(req_packages) > 0) {
  install.packages(req_packages, dependencies = TRUE)
}
sapply(pcakages, require, character.only = TRUE)

# defining working directory
setwd("D:/Marlon/Variables_processing/") # change this to your working directory

# reading data
varaibles_list <- list.files(path = "bio", pattern = ".asc", # vector of variables
                             full.names = TRUE)

variables <- stack(varaibles_list) # stack of variables

# getting data from the variables
variables_values <- na.omit(values(variables))

# sample of 10000 values if more pixels exist
if (dim(variables_values > 10000)) {
  variables_values <- variables_values[sample(rownames(variables_values), 10000), ] 
} 

# correlation matrix calculation
correlation_matrix <- cor(variables_values)

# saving correlation matrix
write.csv(correlation_matrix, "variables_correlation_matrix.csv",
          row.names = FALSE)
