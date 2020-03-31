#########
# Final models and projection, after model calibration
#########

# Description
## The following script helps to create final models and their projections after
## model selection during calibration. To see how model calibration was done, visit
## https://github.com/marlonecobos/ENM_manuals/blob/master/ENM_process/Model_calibration.R

## No data is needed if model calibration was performed as indicated in the link
## above.

## The main processes are performed using the package kuenm from GitHub. To install
## this package see instructions in https://github.com/marlonecobos/kuenm.


# loading needed package
# assuming that you installed kuenm, load it, if not installed see 
# https://github.com/marlonecobos/kuenm for instructions
library(kuenm)


#######################################################################################
# Preparing directory and data
##################

# The data and results necessary to perform the analyses can be obtained following 
# instructions in https://github.com/marlonecobos/kuenm or using the script in
# https://github.com/marlonecobos/ENM_manuals/blob/master/ENM_process/Model_calibration.R

setwd("YOUR/DIRECTORY") # set the working directory


#######################################################################################
# Final models and projections
##################

# Final models without projections
## In this case only models for the area of calibration will be performed.


## check the functions help to understand arguments
help(kuenm_mod)

## preparing arguments 
occ_joint <- "aame_joint.csv"
M_var_dir <- "M_variables"
out_eval <- "Calibration_results"
batch_fin <- "Final_models_code"
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- TRUE
mod_dir <- "Final_Models"
out_format <- "logistic"
project <- FALSE
maxent_path <- "YOUR/DIRECTORY" # where Maxent is
wait1 <- FALSE
run1 <- TRUE
args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
             # "outputgrids=false" which avoids writing grids of replicated models and only writes the 
             # summary of them (e.g., average, median, etc.) when rep.n > 1
             # note that some arguments are fixed in the function and should not be changed

# running final models (no projections)
kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, 
          batch = batch_fin, rep.n = rep_n, rep.type = rep_type, jackknife = jackknife, 
          out.dir = mod_dir, out.format = out_format, project = project,
          maxent.path = maxent_path, args = args, wait = wait1, run = run1)



# Final models with projections
## In this case models will be projected in time and space to the whole world in
## disctinct future scenarions, represented by two RCPs and different GCMs.

## preparing arguments 
occ_joint <- "aame_joint.csv"
M_var_dir <- "M_variables"
out_eval <- "Calibration_results"
batch_fin <- "Final_models_code1"
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- TRUE
mod_dir <- "Final_Models1"
out_format <- "logistic"
project <- TRUE
G_var_dir <- "G_variables"
ext_type <- "all"
write_mess <- FALSE
write_clamp <- FALSE
maxent_path <- "YOUR/DIRECTORY" # where Maxent is
wait1 <- FALSE
run1 <- TRUE
args <- NULL 

# running final models (with projections)
kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, 
          batch = batch_fin, rep.n = rep_n, rep.type = rep_type, jackknife = jackknife, 
          out.dir = mod_dir, out.format = out_format, project = project,
          G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, 
          write.clamp = write_clamp, maxent.path = maxent_path, args = args, 
          wait = wait1, run = run1)