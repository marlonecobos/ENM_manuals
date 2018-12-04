#####
# Creating sets of all potential combinations of a group of variables
#####

# data needed 
## Environmental raster layers, in this example we will use a set of ascii 
## layers located in the folder (bio) in our working directory. Variables must 
## be in ascii format.

## Note: All variables must have the same projection, extent, and resolution.
## We suggest to work with Geographic projections WGS84, with no planar projection.

source("https://raw.githubusercontent.com/marlonecobos/ENM_manuals/master/Variables_processing/kuenm_varcomb.R")

in_folder <- "bio" # name of folder with variables to be combined in distinct sets
out_folder <- "Var_combinations" # name of folder that will contain the sets 
                                 # of combined variables (the function creates this)
n_min <- 2 # minimum number of variables per combination
  
kuenm_varcomb(vars = in_folder, out.vars = out_folder, min.number = n_min)