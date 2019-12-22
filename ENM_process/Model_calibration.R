#########
# Model calibration using kuenm and Maxent
#########

# Description
## The following script helps to perform the process of model calibration using
## Maxent through the R pacakge kuenm. Model calibration consists of a series of
## steps starting in the creation of several models and ending in the selection 
## of parameter settings that produce the best models. To detect which ones are
## the models with best results a robust evaluation process must be performed. In
## this example, models are selected based on statistical significance, omission
## rates based on a predefined Error, and model complexity. The metrics used here,
## are partial ROC, omission rates, and AICc, in that order.

## No data is needed if internet conection is available.

## The main processes are performed the package kuenm from GitHub. To install
## this package see instructions in https://github.com/marlonecobos/kuenm.


# loading needed package
# assuming that you installed kuenm, load it, if not installed see 
# https://github.com/marlonecobos/kuenm for instructions
library(kuenm)


#######################################################################################
# Preparing directory and data
##################

# Data used as an example for testing this package correspond to the turkey tick 
# Amblyomma americanum, a vector of various diseases, including human monocytotropic 
# ehrlichiosis, canine and human granulocytic ehrlichiosis, tularemia, and southern 
# tick-associated rash illness.

# These data are already structured as needed for doing analysis with this package, 
# and can be downloaded (from kuenm example data) and extracted using the code below.

# If you have your own data and they are organized as in the first part of Figure 1, 
# in https://github.com/marlonecobos/kuenm change your directory and avoid the 
# code bellow.

# Change "YOUR/DIRECTORY" by your actual directory.
download.file(url = "https://kuscholarworks.ku.edu/bitstream/handle/1808/26376/ku.enm_example_data.zip?sequence=3&isAllowed=y", 
              destfile = "YOUR/DIRECTORY/ku.enm_example_data.zip", mode = "wb",
              quiet = FALSE) # donwload the zipped example folder in documents

unzip(zipfile = "YOUR/DIRECTORY/ku.enm_example_data.zip",
      exdir = "YOUR/DIRECTORY") # unzip the example folder in documents

unlink("YOUR/DIRECTORY/ku.enm_example_data.zip") # erase zip file

setwd("YOUR/DIRECTORY/ku.enm_example_data/A_americanum") # set the working directory


#######################################################################################
# Complete process of model calibration
##################

# Candidate model creation
## check the functions help to understand arguments
help(kuenm_cal)

## preparing arguments (Change "YOUR/DIRECTORY" by your pertinent directory)
occ_joint <- "aame_joint.csv"
occ_tra <- "aame_train.csv"
M_var_dir <- "M_variables"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
reg_mult <- c(seq(0.1, 1, 0.3), seq(2, 5, 1))
f_clas <- c("lq", "lp", "q", "qp", "lqp")
args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
             # note that some arguments are fixed in the function and should not be changed
maxent_path <- "YOUR/DIRECTORY" # where Maxent is
wait <- FALSE
run <- TRUE

## runing candidate models
kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, 
          batch = batch_cal, out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, 
          args = args, maxent.path = maxent_path, wait = wait, run = run)



# Candidate model evaluation and selection (YOU CAN DO THIS WHILE CANDIDATE MODELS ARE BEING CREATED)
## check the functions help to understand arguments
help(kuenm_ceval)

## preparing arguments (Change "YOUR/DIRECTORY" by your pertinent directory)
occ_test <- "aame_test.csv"
out_eval <- "Calibration_results"
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- TRUE
selection <- "OR_AICc"

## runing candidate models
cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, 
                        occ.test = occ_test, batch = batch_cal, out.eval = out_eval, 
                        threshold = threshold, rand.percent = rand_percent, 
                        iterations = iterations, kept = kept, selection = selection)

# CHECK RESULTS IN WORKING DIRECTORY AND IN THE OBJECT CREATED BEFORE.
# NEXT STEP IS CREATION OF FINAL MODELS
