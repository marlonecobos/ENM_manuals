#####
# Initial cleaning of occurrence data
## Excluding: records with no coordinates, duplicates, records with 
## corrdinates (0, 0), and records with low precision
#####

# data needed 
## (1) Species occurrences; in this example we will use a file named sp_occ.csv. 
## This file must contain the following columns (in that order): 
## ID, Species_name, Longitud, Latitud.

# defining working directory
setwd("C:/Users/Marlon/Documents/R/Data_cleaning") # change this to your working directory

# reading data
occurrences <- read.csv("Cynomys_ludovicianus_Exercise.csv") # occurrences

# Excluding records with no coordinates
occurrences <- occurrences[!is.na(occurrences$long) | !is.na(occurrences$lat), ]

# Excluding duplicates
occurrences$code <-  paste(occurrences$Species, occurrences$long, # concatenating columns of interest
                      occurrences$lat, sep = "_")

occurrences <- occurrences[!duplicated(occurrences$code), 1:4] # erasing duplicates

# Excluding records with (0, 0) coordinates
occurrences <- occurrences[occurrences$long != 0 & occurrences$lat != 0, ]

# Excluding recors with low level of precision (<= 2 decimals)
## samll function to detect precision 
## (from https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r)
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

occurrences <- occurrences[sapply(occurrences$long, decimalplaces) > 2 & # keep only the ones with more than 2 decimals
                             sapply(occurrences$lat, decimalplaces) > 2, ]
