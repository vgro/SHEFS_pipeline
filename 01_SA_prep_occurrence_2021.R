################################################################
# load and prepare occurrence records
# create background/pseudoabsence points
# extract environmental data
#
# author: Vivienne Groner
# latest update: 09.03.2022
###############################################################

# Sourcing and directories
wd<-setwd(my_path)

source(paste0(wd,"05_code_vivienne/pipeline_and_analysis/SA_functions_CBER.R"))
memory.limit(size=500000000)


# Load Libraries
set.seed(2237)
library(CoordinateCleaner)
library(dplyr)
library(here)
library(lwgeom)
library(raster)
library(sp)
library(readr)
library(sf)
library(dismo)
library(rgdal)  


# load reference climate data
climate_path<-paste(wd,"/02_environmental_data/CHELSA/",sep='')

xmin <- 0
xmax <- 45
ymin <- -45
ymax <- 0

ext_occ <- raster::extent(xmin, xmax, ymin, ymax)
min_occ <-  15 # minimum number of occurrence points needed for the species to be included in the models

bio_layers <- list.files(climate_path, pattern = 'tif')
bio_layers1 <- paste(climate_path, bio_layers, sep="/")
bio_layers2 <- lapply(bio_layers1,raster::raster)
env_layers <- raster::stack(bio_layers2)
env_crop <- raster::crop(env_layers, ext_occ)


# load species occurrences from GBIF
GBIF_files <- list.files(paste(wd,'/04_occurrence_records/GBIF_raw/', sep=''),pattern='csv')
foldernames <- gsub('.csv','',GBIF_files)
species_list <- list()


# loop over groups (amphibia, annelida, arthropoda, birds, gastropoda, mammalia, reptilia)
  for (i in 1:length(GBIF_files)){
    
    if(!dir.exists(paste(wd,'/04_occurrence_records/GBIF_raw/',foldernames[[i]],sep=''))){
    dir.create(paste(wd,'/04_occurrence_records/GBIF_raw/',foldernames[[i]],'/',sep=''))
    }
    
    raw1 <- read.csv(paste(wd,'/04_occurrence_records/GBIF_raw/',GBIF_files[[i]], sep=''))[,'species']
    raw1[raw1==""] <- NA
    raw2 <- na.omit(raw1)
    
    species_list[[i]]<-raw2[!duplicated(raw2)]
    write.csv(species_list[[i]],paste(wd,'/04_occurrence_records/GBIF_raw/species_lists/species_list_',GBIF_files[[i]],sep=''))

    unlisted <- gsub(' ','_',species_list[[i]])
    spxy_out <- lapply(X = unlisted, 
                       FUN = gbifData, 
                       ext_sp = ext_occ, 
                       ext_occ = ext_occ, 
                       out_dir = paste(wd, "/04_occurrence_records/GBIF_raw/", 
                                       foldernames[[i]], '/', unlisted,'_',foldernames[[i]],'.csv',sep=''), 
                       min_occ = min_occ)
    
}


# load pre-selected data (species with min occurence at this stage)
GBIF_files<-list.files(paste(wd,'/04_occurrence_records/GBIF_raw/', sep=''),pattern='csv')
foldernames<-gsub('.csv','',GBIF_files)


f=1 # select which group of species
raw_data_dir <- paste(wd,'/04_occurrence_records/GBIF_raw/',foldernames[[f]],sep='')
sp_names <- gsub(".csv", "", list.files(raw_data_dir))
sp_names


# Coordinate Cleaner
if(!dir.exists(paste(wd,"/04_occurrence_records/GBIF_ready/points/cleaned_raw",sep=''))){
  dir.create(paste(wd,"/04_occurrence_records/GBIF_ready/points/cleaned_raw/",sep=''))
}

lapply(X = sp_names, 
       FUN = cc_wrapper, 
       in_dir = raw_data_dir, 
       out_dir = paste(wd,"/04_occurrence_records/GBIF_ready/points/cleaned_raw/",sep=''), 
       min_occ = min_occ)


# Rarefy Points
if(!dir.exists(paste(wd,"/04_occurrence_records/GBIF_ready/points/rarefied",sep=''))){
  dir.create(paste(wd,"/04_occurrence_records/GBIF_ready/points/rarefied/",sep=''))
}

ref_map <- env_crop[[1]]
ref_map[!is.na(ref_map)] <- 0   #ref_map should be full of non-1 values

sp_names <- gsub(".csv", "", list.files(paste(wd,"/04_occurrence_records/GBIF_ready/points/cleaned_raw/",sep=''),pattern='csv'))

lapply(X = sp_names, 
       FUN = rarefyPoints,
       in_dir = paste(wd,"/04_occurrence_records/GBIF_ready/points/cleaned_raw/",sep=''), 
       out_dir = paste(wd,"/04_occurrence_records/GBIF_ready/points/rarefied/",sep=''), 
       ref_map = ref_map, 
       min_occ = min_occ)


# Extract Data for presence points
if(!dir.exists(paste(wd,"/04_occurrence_records/GBIF_ready/environmental/presence/",sep=''))){
  dir.create(paste(wd,"/04_occurrence_records/GBIF_ready/environmental/presence/",sep=''), recursive = TRUE)
}

sp_names <- gsub(".csv", "", list.files(paste(wd,"/04_occurrence_records/GBIF_ready/points/rarefied/",sep=''),pattern='csv'))

lapply(X= sp_names, 
       FUN = ras_extract,
       in_dir = paste(wd,"/04_occurrence_records/GBIF_ready/points/rarefied/",sep=''), 
       out_dir = paste(wd,"/04_occurrence_records/GBIF_ready/environmental/presence/",sep=''), 
       raster_in = env_crop)


# Background Data - create background and pseudoabsence points
if(!dir.exists(paste(wd,"/04_occurrence_records/GBIF_ready/points/background/",sep=''))){
  dir.create(paste(wd,"/04_occurrence_records/GBIF_ready/points/background/",sep=''))
}

if(!dir.exists(paste(wd,"/04_occurrence_records/GBIF_ready/points/pseudoabsence/",sep=''))){
  dir.create(paste(wd,"/04_occurrence_records/GBIF_ready/points/pseudoabsence/",sep=''))
}

ecoreg <-
  sf::st_read(paste(wd,"/02_environmental_data/WWF_Ecoregions/wwf_terr_ecos.shp",sep='')) %>% 
  sf::st_crop(.,ext_occ) %>%  ##cropping to the area of interest
  dplyr::select(OBJECTID, ECO_NAME) ##just selecting out the columns we're interested in
  
lapply(X = sp_names, 
       background_sampler, 
       in_dir= paste(wd,"/04_occurrence_records/GBIF_ready/points/rarefied/",sep=''), 
       out_dir = paste(wd,"/04_occurrence_records/GBIF_ready/points/background/",sep=''), 
       dens_abs = "density", density = 250, type = "background",  polygon = ecoreg)

lapply(X = sp_names, 
       background_sampler, 
       in_dir = paste(wd,"/04_occurrence_records/GBIF_ready/points/rarefied/",sep=''), 
       out_dir = paste(wd,"/04_occurrence_records/GBIF_ready/points/pseudoabsence/",sep=''), 
       dens_abs = "density", density = 250, type = "pseudoabsence", buffer = 25, polygon = ecoreg)


# Extract environmental data for background points 
 if(!dir.exists(paste(wd,"/04_occurrence_records/GBIF_ready/environmental/background/",sep=''))){
   dir.create(paste(wd,"/04_occurrence_records/GBIF_ready/environmental/background/",sep=''))
 }
 
 lapply(X = sp_names, 
        FUN = ras_extract,
        in_dir = paste(wd,"/04_occurrence_records/GBIF_ready/points/background/",sep=''), 
        out_dir = paste(wd,"/04_occurrence_records/GBIF_ready/environmental/background/",sep=''), 
        raster_in = env_crop)


# Extract environmental data for pseudoabsence points
if(!dir.exists(paste(wd,"/04_occurrence_records/GBIF_ready/environmental/pseudoabsence/",sep=''))){
  dir.create(paste(wd,"/04_occurrence_records/GBIF_ready/environmental/pseudoabsence/",sep=''))
}

lapply(X = sp_names, 
       FUN = ras_extract,
       in_dir = paste(wd,"/04_occurrence_records/GBIF_ready/points/pseudoabsence/",sep=''), 
       out_dir = paste(wd,"/04_occurrence_records/GBIF_ready/environmental/pseudoabsence/",sep=''), 
       raster_in = env_crop)


# plot points on map
sp_names <- gsub("*.csv","",list.files(paste(wd,"/04_occurrence_records/GBIF_ready/environmental/presence/",sep=''),
                                       pattern = "*.csv"))

points_in_list <- list.files(paste(wd,"/04_occurrence_records/GBIF_ready/environmental/presence/",sep=''),
                           pattern = "*.csv", full.names = TRUE)
points_in <- lapply(points_in_list,read.csv)

for (i in 1:length(points_in)){
  
  tiff(paste0(paste0(paste(wd,"/04_occurrence_records/GBIF_ready/plots_points/",sep=''), sp_names[[i]], ".png",sep='')))
  plot(ref_map,legend = FALSE, main=paste('occurrence - ', sp_names[[i]]),xlab="longitude", ylab="latitude",col='grey')
  box()
  points(points_in[[i]]$x,points_in[[i]]$y,pch=16)  
  dev.off()
  print(paste(sp_names[[i]]))
  
}

# END