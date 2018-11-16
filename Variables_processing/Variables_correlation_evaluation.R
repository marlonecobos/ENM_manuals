#####
# Evaluation of variables correlation
#####

# data needed 
## Environmental raster layers, in this example we will use a set of ascii 
## layers located in the folder (bio) in our working directory. Our layers are 
## 15 bioclimatic varaibles. These varaibles can be in other formats as well,
## but make sure that you use the appropriate extension when reading them.

## Note: All variables must have the same projection, extent, and resolution.
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

# sample of 10000 values if more pixels exist (optional)
if (dim(variables_values)[1] > 10000) {
  variables_values <- variables_values[sample(1:nrow(variables_values), 10000), ] 
} 

# correlation matrix calculation
correlation_matrix <- cor(variables_values)

# saving correlation matrix
write.csv(correlation_matrix, "variables_correlation_matrix.csv",
          row.names = TRUE)

# detecting correlated varaibles more easily
correlation_matrix1 <- correlation_matrix # making other table with results

max_cor <- 0.8 # maximum value of correlation allowed

for (i in 1:dim(correlation_matrix1)[2]) { #correlated values will turn into 2 for easier detection
  for (j in 1:dim(correlation_matrix1)[1]) {
    correlation_matrix1[j, i] <- ifelse(correlation_matrix1[j, i] > max_cor | correlation_matrix1[j, i] < -max_cor, 
                                        2, correlation_matrix1[j, i])
  }
}

# #checking the table
View(correlation_matrix1) # selection should be done manually, 2 = correlated

# saving correlation matrix
write.csv(correlation_matrix1, "variables_correlation_matrix2.csv",
          row.names = TRUE)

# selecting variables and writing them in a new directory
names(variables) # names

selected_variables <- variables[[c("bio_1", "bio_4", "bio_7", # select only non-correlated ones
                                   "bio_12", "bio_15", "bio_17")]] 

variable_names <- paste("Non_correlated_variables/", names(selected_variables), # names of variables
                        ".asc", sep = "") 

dir.create("Non_correlated_variables") # create the folder for saving these variables

for (i in 1:length(raster::unstack(selected_variables))) { # writing the selected variables
  raster::writeRaster(selected_variables[[i]], filename = variable_names[i], 
                      format = "ascii")
}
