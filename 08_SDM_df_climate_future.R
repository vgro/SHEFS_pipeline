################################################################
# read SDM output and calculate ensemble majority presence/absence
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
library(data.table)

path<-my_path # base directory
path_sdm<-'11_SDM_output/data/' # SDM output
path_LUI<-'02_environmental_data/SSP_landcover/' # land use intensity maps

source(paste(path,'05_code_vivienne/pipeline_and_analysis/SA_functions_CBER_analysis.R',sep=''))

climate_scen_name<-c('rcp45_ACCESS1-0','rcp45_CNRM-CM5','rcp45_GFDL-ESM2G','rcp45_HadGEM2-CC','rcp45_MPI-ESM-MR',
                     'rcp85_ACCESS1-0','rcp85_CNRM-CM5','rcp85_GFDL-ESM2G','rcp85_HadGEM2-CC','rcp85_MPI-ESM-MR')


# South Africa extent and shape
extent<-c(16.4, 34, -35, -21.5)
mask_SA_in<-readOGR(paste(path,'02_environmental_data/Igismap/SouthAfrica_Bound.shp',sep=''))
mask_SA<-crop(mask_SA_in,extent)


# get species names from directories or list
dirs<-list.dirs(paste(path,path_sdm,sep=''),full.names = FALSE, recursive = FALSE)
species_names<-dirs
species_names


# loop over climate scenarios
for (fut in 1:length(climate_scen_name)){
  
  #filename for data frame climatic habitat suitability
  outfile<-paste(path,'11_SDM_output/SA_data_', climate_scen_name[[fut]],'.RDS',sep='') 
  print(climate_scen_name[[fut]])
  
  # create empty data.frame for presences with species name, lat lon, scenario
  df_climate<-as.data.frame(matrix(ncol=4))
  colnames(df_climate)<-c('species','x','y',climate_scen_name[[fut]])

  ### loop over species
  for (s in 1:length(species_names)){
    sp_names<-species_names[[s]]
    base_dir_out<-paste(path,path_sdm,sp_names,'/',sep='')
  
    file<-paste(paste(base_dir_out,'future/',climate_scen_name[[fut]],'/predictions/ensemble/majority_pa/',sp_names,'_ensemble.tif',sep=''))
    
    if (file.exists(file)){
      
      # climate only
      HS_maps_future<-raster(file)
      HS_maps_future_cr<-crop(HS_maps_future,extent)
      HS_maps_future_c<-mask(HS_maps_future_cr,mask_SA)
    
      if (cellStats(HS_maps_future_c,'sum')==0){
        print(paste(sp_names, 'not present in scenario', climate_scen_name[[fut]],sep=' '))
      }    
    
      else {
        xyf<-as.data.frame(HS_maps_future_c,xy=TRUE) # raster to data frame 
        pres_cellsf<-as.data.frame(c(sp_names,subset(xyf,xyf[,3]==1)))
        colnames(pres_cellsf)<-c('species','x','y',paste(climate_scen_name[[fut]]))
        l<-list(df_climate,pres_cellsf)
        df_climate<-data.table::rbindlist(l,use.names=TRUE)
        #print(paste(sp_names, climate_scen_name[[fut]], 'complete',sep=' '))
      }
    }
    
    else { # if file does not exist
      
    #### Get evaluation and raw output to calculate majority pa
      eval_files <-
        list.files(
          paste0(base_dir_out,"evaluation/"),
          full.names = TRUE,
          recursive = TRUE,
          pattern = "*eval")
      
      evals_out <- lapply(eval_files, get_eval, threshold = "spec_sens")
      eval_df <- do.call(rbind, evals_out)
      eval_df$sp_name <- as.character(eval_df$sp_name)
      
      predsf <- list.files(paste0(base_dir_out, "future/",climate_scen_name[[fut]],"/predictions"), full.names = TRUE, recursive = TRUE)
      predsf <- predsf[!grepl("/ensemble/", predsf)]
      
      sp_outf <- lapply(sp_names, get_spf)
      eval_namesf <- do.call(rbind, sp_outf)  
      
      # climate only
      HS_maps_future<-ensemble_model_f(sp_names, eval_df = eval_df, preds = predsf, method = "majority_pa")
      HS_maps_future_cr<-crop(HS_maps_future,extent)
      HS_maps_future_c<-mask(HS_maps_future_cr,mask_SA)
      
      if (cellStats(HS_maps_future_c,'sum')==0){
        print(paste(sp_names, 'not present in scenario', climate_scen_name[[fut]],sep=' '))
      }    
      
      else {
        xyf<-as.data.frame(HS_maps_future_c,xy=TRUE) # raster to data frame 
        pres_cellsf<-as.data.frame(c(sp_names,subset(xyf,xyf[,3]==1)))
        colnames(pres_cellsf)<-c('species','x','y',paste(climate_scen_name[[fut]]))
        l<-list(df_climate,pres_cellsf)
        df_climate<-data.table::rbindlist(l,use.names=TRUE)
        #print(paste(sp_names, climate_scen_name[[fut]], 'complete',sep=' '))
        }
     }
  }
  
  df_climate<-df_climate[rowSums(is.na(df_climate[,2:4])) < 3,]
  saveRDS(df_climate,outfile)
}


# END

