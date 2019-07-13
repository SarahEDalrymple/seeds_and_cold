# Freezing tolerance and seed structure

# Sarah Dalrymple
# May 2019

# An analysis to determine the coldest temperature at which species of plants can survive in habitat

# Workflow

# 1. download climate data
# 2. download species occurrence data
# 3. determine the location at which the focal species experiences the lowest temperature
# 4. derive climate data for that point

##############################################################################################

### Step 0: preparing R for analysis
#########################################

# load libraries needed for analysis
# if not yet installed, R will return an error message in the console,
# use install.packages() inserting the name of the package in quotation marks in the brackets,
# and then try loading the libraries again

library(rgbif)
library(raster)
library(rgdal)

#######################################################################################

# 1.download climate data
#########################

# from http://worldclim.org/version2 download the full set of bioclim variables at 2.5arcminute resolution,
# this is equivalent to 4.5km at the equator

#load the .tif files as rasters

#BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#BIO3 = Isothermality (BIO2/BIO7) (* 100)
#BIO4 = Temperature Seasonality (standard deviation *100)
#BIO6 = Min Temperature of Coldest Month
#BIO7 = Temperature Annual Range (BIO5-BIO6)

file.choose()

diurnal <- raster("C:\\Users\\sarah\\Dropbox\\Worldclim downloads\\wc2.0_30s_bio\\wc2.0_bio_30s_02.tif")
isotherm <- raster("C:\\Users\\sarah\\Dropbox\\Worldclim downloads\\wc2.0_30s_bio\\wc2.0_bio_30s_03.tif")
t_season <- raster("C:\\Users\\sarah\\Dropbox\\Worldclim downloads\\wc2.0_30s_bio\\wc2.0_bio_30s_04.tif")
t_min <- raster("C:\\Users\\sarah\\Dropbox\\Worldclim downloads\\wc2.0_30s_bio\\wc2.0_bio_30s_06.tif")
t_ann_range <- raster("C:\\Users\\sarah\\Dropbox\\Worldclim downloads\\wc2.0_30s_bio\\wc2.0_bio_30s_07.tif")


climate <- stack(diurnal,isotherm, t_season, t_min, t_ann_range)
plot(climate)
#######################################################################################

# 2. download species occurrence data
######################################


wd <- setwd ("C:\\Users\\sarah\\Dropbox\\Seeds\\Alpine seeds\\Endo_vs_nonendo\\Data_extraction_30stest")

species <- read.table("C:\\Users\\sarah\\Dropbox\\Seeds\\Alpine seeds\\Endo_vs_nonendo\\Species15.txt", header = TRUE)

head(species)

speciesList <- as.list(species$name)


for(i in speciesList) {searchSpecies <- as.character(i)

searchSpeciesAccepted <- name_suggest(q= i, rank='species')

#if(nrow(searchSpeciesAccepted) == 0) {print(paste("no match", i))
# setwd ("C:\\Users\\sarah\\Dropbox\\Seeds\\Alpine seeds\\Endo_vs_nonendo\\Data_extraction_30s")

#  no_match <- data.frame("species" = i, 
#                        "key" = "NA", 
#                        "Name" = "NA", 
#                        "Country" = "NA")
#  class(no_match)
# write.table(no_match,file="all_species_dataextraction1.csv", 
# append=TRUE, sep= ",", row.names = FALSE, col.names=FALSE)}

key <- name_backbone(name= searchSpecies)$speciesKey

no_occ <- occ_count(georeferenced = TRUE, taxonKey = key)

if(no_occ > 30000) {print(paste("return limit exceeded", i))
  setwd ("C:\\Users\\sarah\\Dropbox\\Seeds\\Alpine seeds\\Endo_vs_nonendo\\Data_extraction_30stest")
  
  return_limit_exceeded <- data.frame("species" = i, 
                                      "key" = key, 
                                      "Name" = searchSpecies, 
                                      "Country" = "NA")
  class(return_limit_exceeded)
  write.table(return_limit_exceeded,file="all_species_dataextraction.csv", 
              append=TRUE, sep= ",", row.names = FALSE, col.names=FALSE)
} else {
  
  searchSpeciesUnderscore <- sub(" ", "_", searchSpecies)
  searchSpeciesUnderscore
  newdir<- dir.create(searchSpeciesUnderscore)
  setwd(file.path(searchSpeciesUnderscore))
  
  
  occ <- occ_search(scientificName =  paste(searchSpecies),
                    fields = c('key', 'scientificName','country','decimalLatitude','decimalLongitude'),
                    hasCoordinate=T ,
                    limit = 50000,
                    return = 'data')
  
  #remove blank spaces from species names.
  occ$scientificName <- sub(" ", ".", occ$scientificName)
  
  
  #####################################################################################
  
  # 3. determine the location experiencing lowest temperature
  
  
  coordinates(occ)= ~ decimalLongitude+ decimalLatitude
  #  plot(occ)
  
  # if you want to view points against a base map use the following code substituting the file path for your own
  # however, it can take some time to draw the map so this stage can be skipped without affecting data extraction
  #global_map <-readOGR("C:\\Users\\sarah\\Dropbox\\GIS_inc_DIVA\\Global map files\\countries.shp")
  #plot(global_map)
  #points(occ, pch = 16, col = "red")
  
  # Extract raster value by points
  
  min_temp_occ = extract(t_min, occ)
  
  # Combine raster values with point and save as a CSV file.
  
  combinePointValue=cbind(occ,min_temp_occ)
  
  combinePointValue <- data.frame(combinePointValue)
  colnames(combinePointValue) <- c("key", "name", "country", "t_min_value", 
                                   "decimalLongitude", "decimalLatitude","optional")
  
  # write.table below only needed if data needs to be stored
  write.table(combinePointValue,file="t_min_AllOccurrences.csv", 
              append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)
  
  min_temp_loc <- combinePointValue[order(combinePointValue$t_min_value),]
  min_temp_loc <- min_temp_loc[1,]
  
  #######################################################################################
  
  # 4. extract other climate variables for lowest temp location
  
  coordinates(min_temp_loc)= ~ decimalLongitude+ decimalLatitude
  
  
  # Extract raster value by points
  
  rasValues=extract(climate, min_temp_loc)
  
  # Combine raster values with point and save as a CSV file.
  combinePointValue=cbind(min_temp_loc,rasValues)
  write.table(combinePointValue,file="min_temp_loc combinedPointValue.csv", 
              append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE)
  
  ########################################################################################
  
  # 5. append climate data to full dataset
  
  setwd ("C:\\Users\\sarah\\Dropbox\\Seeds\\Alpine seeds\\Endo_vs_nonendo\\Data_extraction_30stest")
  
  
  write.table(combinePointValue,file="all_species_dataextraction.csv", 
              append=TRUE, sep= ",", row.names = i, col.names=FALSE)
  
  #########################################################################################
  }
  }
