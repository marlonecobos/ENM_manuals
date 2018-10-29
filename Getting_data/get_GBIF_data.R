#########
# Download species information from GBIF using rgbif
#########

# Description
## The following script helps to download species occurrence data from the GBIF
## database in three distinct ways: (1) for a single species; (2) for multiple
## species not necessaly related to each other, and (3) for all species of one or
## multiple genus.

## In addition, two first steps of data cleaning are performed when downloading after
## downoading the data, these are: erasing records without georeference and elminating
## duplicates. If these steps are not needed, add a # beofer the lines that are 
## idicated in the script


# loading needed package
if(!require(rgbif)){
  install.packages("rgbif")
  library(rgbif)
}


#######################################################################################
# Single species #
##################

# defining working directory
## project forlder
setwd("Z:/Marlon_E_Cobos/ENM_project") # Your folder

# getting the data from GBIF
sp <- "Dasypus kappleri" # species name

sps <- name_lookup(query = sp, rank = "species", return = "data") # information about the species

keys <- sps$key # all keys returned
counts <- vector() # object to save info on number of records per key

for (j in 1:length(keys)) { # testing if keys return records
  counts[j] <- occ_count(taxonKey = keys[j], georeferenced = TRUE) 
}

if (sum(counts) == 0) { # if no info, tell the species
  cat("species", sp, "has no goereferenced data\n")
  
}else { # if it has info, use the key with more records, which is the most useful
  if (length(keys) == 1) { # if it is only one key
    key <- keys # detecting species key 

  }else { # if its more than one key
    keysco <- cbind(keys, counts)
    keysco <- keysco[order(keysco[, 2]), ]
    key <- keysco[dim(keysco)[1], 1] # detecting species key that return information
  }
  
  occ <- occ_search(taxonKey = key, return = "data", limit = 10000) # getting the data from GBIF
  occ_g <- occ
  occ_g <- occ_g[, c(1, 2, 4, 3, 5:dim(occ_g)[2])] # reordering longitude and latitude
  
  # keeping only unique georeferenced records. IF NO FILTERING IS NEEDED, PUT A # IN FRONT OF THE NEXT 3 LINES
  occ_g <- occ_g[!is.na(occ_g$decimalLatitude) & !is.na(occ_g$decimalLongitude), ] # excluding no georeferences
  occ_g <- occ_g[!duplicated(paste(occ_g$name, occ_g$decimalLatitude, # excluding duplicates
                                   occ_g$decimalLongitude, sep = "_")), ]
  
  # writting file
  file_name <- paste(gsub(" ", "_", sp), "csv", sep = ".") # csv file name per each species
  write.csv(occ_g, file_name, row.names = FALSE) # writing inside each genus folder
}



#######################################################################################
# Multiple species #
####################

# defining working directory
## project forlder
setwd("Z:/Marlon_E_Cobos/ENM_project") # Your folder

spvector <- c("Dasypus kappleri", "Panthera onca", "Dasyprocta punctata") # binomial names

## Getting info species by species 
occ_count <- list() # object to save info on number of georeferenced records per species 

for (i in 1:length(spvector)) {
  sps <- name_lookup(query = spvector[i], rank = "species", 
                     return = "data", limit = 100) # information about the species
  
  keys <- sps$key # all keys returned
  counts <- vector() # object to save info on number of records per key
  
  for (j in 1:length(keys)) { # testing if keys return records
    counts[j] <- occ_count(taxonKey = keys[j], georeferenced = TRUE) 
  }
  
  if (sum(counts) == 0) { # if no info, tell the species
    cat("species", spvector[i], "has no goereferenced data\n")
    
  }else { # if it has info, use the key with more records, which is the most useful
    if (length(keys) == 1) { # if it is only one key
      key <- keys # detecting species key 
      occ_count[[i]] <- cbind(spvector[i], counts) # count how many records
      
    }else { # if its more than one key
      keysco <- cbind(keys, counts)
      keysco <- keysco[order(keysco[, 2]), ]
      key <- keysco[dim(keysco)[1], 1] # detecting species key that return information
      occ_count[[i]] <- keysco[dim(keysco)[1], ]# count how many records
    }
    
    occ <- occ_search(taxonKey = key, return = "data", limit = 10000) # getting the data from GBIF
    occ_g <- occ
    occ_g <- occ_g[, c(1, 2, 4, 3, 5:dim(occ_g)[2])] # reordering longitude and latitude
    
    # keeping only unique georeferenced records. IF NO FILTERING IS NEEDED, PUT A # IN FRONT OF THE NEXT 3 LINES
    occ_g <- occ_g[!is.na(occ_g$decimalLatitude) & !is.na(occ_g$decimalLongitude), ] # excluding no georeferences
    occ_g <- occ_g[!duplicated(paste(occ_g$name, occ_g$decimalLatitude, # excluding duplicates
                                     occ_g$decimalLongitude, sep = "_")), ]
    
    # writting file
    file_name <- paste(gsub(" ", "_", spvector[i]), "csv", sep = ".") # csv file name per each species
    write.csv(occ_g, file_name, row.names = FALSE) # writing inside each genus folder
    
    cat(i, "of", length(spvector), "species\n") # counting species per genus 
  }
}

occ_count[sapply(occ_count, is.null)] <- NULL # excluding countings with no data
genus_data <- do.call(rbind, occ_count) # making the list of countings a table
genus_data <- data.frame(genus_data[, 1], as.numeric(genus_data[, 2])) # making countings numeric
names(genus_data) <- c("Species", "N_records") # naming columns

# writing the table
file_nam <- "Species_record_count.csv" # csv file name for all species
write.csv(genus_data, file_nam, row.names = FALSE) # writing inside each genus folder



#######################################################################################
# All species of one or more genus #
####################################

# defining working directory
## project forlder
setwd("Z:/Marlon_E_Cobos/ENM_project") # your folder

# defining general variables
genus <- c("Erythranthe", "Mimulus") # use only one genus, or add more if you need

for (h in 1:length(genus)) {
  ## genus folder
  dir.create(genus[h])
  infolder <- paste(getwd(), genus[h], sep = "/")
  
  # genus Erythranthe
  ## all species list
  species <- name_lookup(query = genus[h], rank = "species", 
                         return = "data", limit = 1000) # information about the species
  
  ## working to get unique binomial names
  species_vec <- species$scientificName # vector of scientific names
  
  pattern <- paste(genus[h], "\\S*") # pattern to look for
  species_vect <- gregexpr(pattern, species_vec)
  species_vecto <- regmatches(species_vec, species_vect)
  species_vector <- unlist(species_vecto) # vector of all binomial names found
  
  spvector <- unique(species_vector) # unique binomial names
  
  ## Getting info species by species 
  occ_count <- list() # object to save info on number of georeferenced records per species 
  
  for (i in 1:length(spvector)) {
    sps <- name_lookup(query = spvector[i], rank = "species", 
                       return = "data", limit = 100) # information about the species
    
    keys <- sps$key # all keys returned
    counts <- vector() # object to save info on number of records per key
    
    for (j in 1:length(keys)) { # testing if keys return records
      counts[j] <- occ_count(taxonKey = keys[j], georeferenced = TRUE) 
    }
    
    if (sum(counts) == 0) { # if no info, tell the species
      cat("species", spvector[i], "has no goereferenced data\n")
      
    }else { # if it has info, use the key with more records, which is the most useful
      if (length(keys) == 1) { # if it is only one key
        key <- keys # detecting species key 
        occ_count[[i]] <- cbind(spvector[i], counts) # count how many records
        
      }else { # if its more than one key
        keysco <- cbind(keys, counts)
        keysco <- keysco[order(keysco[, 2]), ]
        key <- keysco[dim(keysco)[1], 1] # detecting species key that return information
        occ_count[[i]] <- keysco[dim(keysco)[1], ]# count how many records
      }
      
      occ <- occ_search(taxonKey = key, return = "data", limit = 10000) # getting the data from GBIF
      occ_g <- occ
      occ_g <- occ_g[, c(1, 2, 4, 3, 5:dim(occ_g)[2])] # reordering longitude and latitude
      
      # keeping only unique georeferenced records. IF NO FILTERING IS NEEDED, PUT A # IN FRONT OF THE NEXT 3 LINES
      occ_g <- occ_g[!is.na(occ_g$decimalLatitude) & !is.na(occ_g$decimalLongitude), ] # excluding no georeferences
      occ_g <- occ_g[!duplicated(paste(occ_g$name, occ_g$decimalLatitude, # excluding duplicates
                                       occ_g$decimalLongitude, sep = "_")), ]
      
      # writting file
      file_name <- paste(gsub(" ", "_", spvector[i]), "csv", sep = ".") # csv file name per each species
      write.csv(occ_g, paste(infolder, file_name, sep = "/"), row.names = FALSE) # writing inside each genus folder
      
      cat("    ", i, "of", length(spvector), "species\n") # counting species per genus 
    }
  }
  
  occ_count[sapply(occ_count, is.null)] <- NULL # excluding countings with no data
  genus_data <- do.call(rbind, occ_count) # making the list of countings a table
  genus_data <- data.frame(genus_data[, 1], as.numeric(genus_data[, 2])) # making countings numeric
  names(genus_data) <- c("Species", "N_records") # naming columns
  
  # writing the table
  file_name1 <- paste(genus[h], "record_count.csv", sep = "_") # csv file name per each genus
  write.csv(genus_data, file_name1, row.names = FALSE) # writing inside each genus folder
  
  cat(h, "of", length(genus), "genus\n") # counting genus ready
}