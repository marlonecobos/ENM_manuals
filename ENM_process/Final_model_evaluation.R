#########
# Final model evaluation
#########

# Description
## The following script helps to evaluate final models when independent data is 
## available. This is done after model calibration, selection, and final model
## creation. To see how model calibration and final models were done, visit
## https://github.com/marlonecobos/ENM_manuals/blob/master/ENM_process/Model_calibration.R
## and https://github.com/marlonecobos/ENM_manuals/blob/master/ENM_process/Final_models.R

## No data is needed if model calibration and final models were performed as 
## indicated in the links above.

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
# Final model evaluation
##################

# Final models should be evaluated using independent occurrence data (i.e., data 
# that have not been used in the calibration process that usually come from 
# different sources). The kuenm_feval function evaluates final models based on 
# statistical significance (partial ROC) and omission rate (E). 

## functions help
help(kuenm_feval)

## preparing arguments
mod_dir <- "Final_Models"
occ_joint <- "aame_joint.csv"
occ_ind <- "aame_ind.csv"
replicates <- TRUE
out_feval <- "Final_Models_evaluation"
threshold <- 5
rand_percent <- 50
iterations <- 500

## running final model evaluation
fin_eval <- kuenm_feval(path = mod_dir, occ.joint = occ_joint, occ.ind = occ_ind, 
                        replicates = replicates, out.eval = out_feval, 
                        threshold = threshold, rand.percent = rand_percent,
                        iterations = iterations)