################################################################
# read SDM output and calculate ensemble majority presence/absence
# for current climate
#
# author: Vivienne Groner
# latest update: 09.03.2022
###############################################################

# set environment and load libraries
library(raster)
library(sp)
library(rgdal)
library(dplyr)
library(tidyverse)

path<-my_path # base directory
path_sdm<-'11_SDM_output/data/' # SDM output
path_LUI<-'02_environmental_data/SSP_landcover/' # land use intensity maps
outfile<-paste(path,'11_SDM_output/SA_data_climate_current.RDS',sep='') # data frame climatic habitat suitability

source(paste(path,'05_code_vivienne/pipeline_and_analysis/SA_functions_CBER_analysis.R',sep=''))


# South Africa extent and shape
extent<-c(16.4, 34, -35, -21.5)
mask_SA_in<-readOGR(paste(path,'02_environmental_data/Igismap/SouthAfrica_Bound.shp',sep=''))
mask_SA<-crop(mask_SA_in,extent)
#plot(mask_SA)


# get species names from directories or list
dirs<-list.dirs(paste(path,path_sdm,sep=''),full.names = FALSE, recursive = FALSE)
species_names<-dirs


# create empty data.frame for presences with species name, lat lon
df_climate<-as.data.frame(matrix(ncol=4))
colnames(df_climate)<-c('species','x','y','current')
print(Sys.time())


### loop over species
for (s in 1:length(species_names)){
  sp_names<-species_names[[s]]
  base_dir_out<-paste(path,path_sdm,sp_names,'/',sep='')
  
# load and crop habitat suitability map
  file<-paste(base_dir_out,'predictions/ensemble/majority_pa/',sp_names,'_ensemble.tif',sep='')
  if (file.exists(file)){
  HS_maps_pres<-raster(file)
  HS_maps_pres_cr<-crop(HS_maps_pres,extent)
  HS_maps_pres_c<-mask(HS_maps_pres_cr, mask_SA)
  
  if (cellStats(HS_maps_pres_c,'sum')==0){
    print(paste(sp_names, 'not present in baseline',sep=' '))
  }  
  
  else {
    xy<-as.data.frame(HS_maps_pres_c,xy=TRUE) # raster to data frame
    pres_cells<-as.data.frame(c(sp_names,subset(xy,xy[,3]==1)))
    colnames(pres_cells)<-c('species','x','y','current')
    l<-list(df_climate,pres_cells)
    df_climate<-data.table::rbindlist(l,use.names=TRUE)
    #print(paste(sp_names, 'present complete',sep=' '))
  }
  }
  else{
    print(paste(sp_names, ' not available',sep = ''))
  }
}

head(df_climate)
df_climate<-df_climate[rowSums(is.na(df_climate[,2:4])) < 3,]
saveRDS(df_climate,outfile)

# END

