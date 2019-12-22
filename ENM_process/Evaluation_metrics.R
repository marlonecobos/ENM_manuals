#########
# Examples of evaluation metrics applied independently
#########

# Description
## The following script helps to run examples of evaluation metrics for ecological
## niche models. Partial ROC and omission rates can be applied to different model
## outputs. AICc can be applied to models from Maxent or the ones in which the
## number of parameters is known, but users have to consider that such metric is
## applied on raster predictions and not in the data used to fit the model.

## No data is needed if the kuenm package is installed.

## The main processes are performed the package kuenm from GitHub. To install
## this package see instructions in https://github.com/marlonecobos/kuenm.


# loading needed package
# assuming that you installed kuenm, load it, if not installed see 
# https://github.com/marlonecobos/kuenm for instructions
library(kuenm)


#######################################################################################
# Preparing directory and data
##################

# defining working directory
## project forlder
setwd("Z:/Marlon_E_Cobos/ENM_project") # Your folder

# getting example data
## train occurrences
octr <- read.csv(list.files(system.file("extdata", package = "kuenm"),
                            pattern = "sp_train.csv", full.names = TRUE))

## test occurrences
occ <- read.csv(list.files(system.file("extdata", package = "kuenm"),
                           pattern = "sp_test.csv", full.names = TRUE))

## all occurrences
ocj <- read.csv(list.files(system.file("extdata", package = "kuenm"),
                           pattern = "sp_joint.csv", full.names = TRUE))

## model to tested partial ROC and omission rates
model <- raster::raster(list.files(system.file("extdata", package = "kuenm"),
                                   pattern = "sp_model.tif", full.names = TRUE))

## model to tested partial ROC and omission rates
modelAICc <- raster::raster(list.files(system.file("extdata", package = "kuenm"),
                                       pattern = "sp_model_joint.tif", full.names = TRUE))

## lambdas file from maxent model
lbds <- readLines(list.files(system.file("extdata", package = "kuenm"), # lambdas file
                             pattern = "lambdas_model_joint.lambdas", full.names = TRUE))



#######################################################################################
# Examples of aplication
##################

# partial ROC
## functions help
help(kuenm_proc)

## arguments
thres <- 5 # percentage of error
rand_perc <- 50 # percentage of data for bootstrapping
iterac <- 500 # number of iterations to get p value

## running
p_roc <- kuenm_proc(occ.test = occ, model = model, threshold = thres,
                    rand.percent = rand_perc, iterations = iterac)



# omission rates
## functions help
help(kuenm_omrat)

## arguments
thres <- 5 # percentage of error

## running
om_rate <- kuenm_omrat(model, threshold = thres, occ.tra = octr, occ.test = occ)



# AICc
## functions help
help(kuenm_aicc)

## arguments
npar <- n.par(lbds) # number of parameters in this model

## running
aicc <- kuenm_aicc(occ = ocj, model = modelAICc, npar = npar)