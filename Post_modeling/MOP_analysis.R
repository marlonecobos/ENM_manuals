#########
# Mobility-oriented parity (MOP) analysis
#########

# Description
## The following script helps to perform the MOP analysis as described in 
## Owens et al. (2013; https://www.sciencedirect.com/science/article/pii/S0304380013002159). 
## To see how the data needed for performing this analysis is organized and gathered 
## check https://github.com/marlonecobos/ENM_manuals/blob/master/ENM_process/Final_models.R

## No data is needed if they were obtained as indicated in the previous link.

## The main processes are performed using the package kuenm from GitHub. To install
## this package see instructions in https://github.com/marlonecobos/kuenm.


# loading needed package
# assuming that you installed kuenm, load it, if not installed see 
# https://github.com/marlonecobos/kuenm for instructions
library(kuenm)
library(raster)


#######################################################################################
# Preparing directory and data
##################

# The data necessary to perform the analyses can be obtained following instructions
# in https://github.com/marlonecobos/kuenm or using the script in
# https://github.com/marlonecobos/ENM_manuals/blob/master/ENM_process/Model_calibration.R

setwd("YOUR/DIRECTORY") # set the working directory


#######################################################################################
# MOP for a single scenario of projection 
##################

## check the functions help to understand arguments
help(kuenm_mod)

## preparing arguments
m_vars <- list.files(path = "M_variables", pattern = ".asc$", full.names = TRUE)
g_vars <- list.files(path = "G_variables/Set3/ccsm4_4.5", pattern = ".asc$", 
                     full.names = TRUE)

m_stack <- stack(m_vars)
g_stack <- stack(g_vars)
percent <- 5


# running final models (no projections)
mop_results <- kuenm_mop(M.stack = m_stack, G.stack = g_stack, percent = percent)



#######################################################################################
# MOP for multiple scenarios of projection 
##################

## preparing arguments 
M_var_dir <- "M_variables"
G_var_dir <- "G_variables"
sets_var <- "Set3" 
out_mop <- "MOP_results"
percent <- 5

# running final models (with projections)
kuenm_mmop(G.var.dir = G_var_dir, M.var.dir = M_var_dir, sets.var = sets_var, 
           out.mop = out_mop, percent = percent)