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
suppressWarnings({
  if(!require(spocc)){
    install.packages("spocc")
    library(spocc)
  }
  #if(!require(scrubr)){
  #  install.packages("scrubr")
  #  library(scrubr)
  #}
  if(!require(rgbif)){
    install.packages("rgbif")
    library(rgbif)
  }
})


#######################################################################################
# Single species #
##################

# defining working directory
## project forlder
setwd("Z:/Marlon_E_Cobos/ENM_project") # Your folder

# getting the data from GBIF
sp <- "Dasypus kappleri" # species name

occ <- occ(query = sp, from = "gbif", limit = 1000) # getting data
occ <- fixnames(occ, how = "query")$gbif$data[[1]] # fix_names

# keeping only unique georeferenced records.
occ_g <- occ[!is.na(occ$longitude) & !is.na(occ$latitude), ] # excluding no georeferences
occ_g <- occ_g[!duplicated(paste(occ_g$name, occ_g$longitude, # excluding duplicates
                                 occ_g$latitude, sep = "_")), ]
occ_g <- occg[, c("name", "longitude", "latitude")] # only these three columns

# writting files
file_name <- paste0(gsub(" ", "_", sp), "_gbif.csv") # csv file name 
write.csv(occ, file_name, row.names = FALSE) # writing

file_name <- paste0(gsub(" ", "_", sp), "_georef.csv") # csv file name 
write.csv(occ_g, file_name, row.names = FALSE) # writing




#######################################################################################
# Multiple species #
####################

# defining working directory
## project forlder
setwd("Z:/Marlon_E_Cobos/ENM_project") # Your folder

spvector <- c("Dasypus kappleri", "Panthera onca", "Dasyprocta punctata") # binomial names

## Getting info species by species 
occ_count <- lapply(1:length(spvector), function(i) {
  occ <- occ(query = spvector[i], from = "gbif", limit = 1000) # getting the data from GBIF
  occ <- fixnames(occ, how = "query")$gbif$data[[1]] # fix_names
  
  # keeping only unique georeferenced records.
  occ_g <- occ[!is.na(occ$longitude) & !is.na(occ$latitude), ] # excluding no georeferences
  occ_g <- occ_g[!duplicated(paste(occ_g$name, occ_g$longitude, # excluding duplicates
                                   occ_g$latitude, sep = "_")), ]
  occ_g <- occg[, c("name", "longitude", "latitude")] # only these three columns
  
  # writting file
  file_name <- paste0(gsub(" ", "_", spvector[i]), "_gbif.csv") # csv file name 
  write.csv(occ, file_name, row.names = FALSE) # writing
  
  file_name <- paste0(gsub(" ", "_", spvector[i]), "_georef.csv") # csv file name 
  write.csv(occ_g, file_name, row.names = FALSE) # writing
  
  cat(i, "of", length(spvector), "species\n") # counting species per genus 
  return(c(Species = spvector[i], counts = nrow(occ_g))) # returning n of records per species 
}) 

count_data <- do.call(rbind, occ_count) # making the list of countings a table
colnames(count_data) <- c("Species", "N_records") # naming columns

# writing the table
file_nam <- "Species_record_count.csv" # csv file name for all species
write.csv(count_data, file_nam, row.names = FALSE) # writing 



#######################################################################################
# All species of one or more genus (STILL UNDER CONSTRUCTION)
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
  occ_count <- lapply(1:length(spvector), function(i) {
    occ <- occ(query = spvector[i], from = "gbif", limit = 1000) # getting the data from GBIF
    occ <- fixnames(occ, how = "query")$gbif$data[[1]] # fix_names
    
    # keeping only unique georeferenced records.
    occ_g <- occ[!is.na(occ$longitude) & !is.na(occ$latitude), ] # excluding no georeferences
    occ_g <- occ_g[!duplicated(paste(occ_g$name, occ_g$longitude, # excluding duplicates
                                     occ_g$latitude, sep = "_")), ]
    occ_g <- occg[, c("name", "longitude", "latitude")] # only these three columns
    
    # writting file
    file_name <- paste0(gsub(" ", "_", spvector[i]), "_gbif.csv") # csv file name 
    write.csv(occ, paste(infolder, file_name, sep = "/"), row.names = FALSE) # writing
    
    file_name <- paste0(gsub(" ", "_", spvector[i]), "_georef.csv") # csv file name 
    write.csv(occ_g, paste(infolder, file_name, sep = "/"), row.names = FALSE) # writing
    
    cat(i, "of", length(spvector), "species\n") # counting species per genus 
    return(c(Species = spvector[i], counts = nrow(occ_g))) # returning n of records per species 
  })
  
  genus_data <- do.call(rbind, occ_count) # making the list of countings a table
  colnames(genus_data) <- c("Species", "N_records") # naming columns
  
  # writing the table
  file_name1 <- paste0(genus[h], "_record_count.csv") # csv file name per each genus
  write.csv(genus_data, file_name1, row.names = FALSE) # writing inside each genus folder
  
  cat(h, "of", length(genus), "genus\n") # counting genus ready
}