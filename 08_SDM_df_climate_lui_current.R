################################################################
# read SDM output and calculate ensemble majority presence/absence
# for current climate and land use and intensity effects
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

source(paste(path,'05_code_vivienne/pipeline_and_analysis/SA_functions_CBER_analysis.R',sep=''))


# South Africa extent and shape
extent<-c(16.4, 34, -35, -21.5)
mask_SA_in<-readOGR(paste(path,'02_environmental_data/Igismap/SouthAfrica_Bound.shp',sep=''))
mask_SA<-crop(mask_SA_in,extent)
plot(mask_SA)


# get species names from directories or list
dirs<-list.dirs(paste(path,path_sdm,sep=''),full.names = FALSE, recursive = FALSE)
species_names<-dirs


# land use and intensity
crop_scen_names<-c('intense')#,'light')
lui_scen_names<-c('X2020_ssp2')

for( lui in 1:length(lui_scen_names)){
  infiles_LUI<-list.files(paste(path,path_LUI,sep=''),pattern=lui_scen_names[[lui]]) # adjust here
  infiles_LUI_r<-lapply(paste(path,path_LUI,infiles_LUI,sep=''),raster)
  infiles_LUI_c<-lapply(infiles_LUI_r,crop,extent)
  files_LUI<-lapply(infiles_LUI_c,mask,mask_SA)
  files_LUI<-stack(files_LUI)

  # loop over different crop intensities
  for (crop in 1:length(crop_scen_names)){
    df_climate_lui<-as.data.frame(matrix(ncol=4))
    colnames(df_climate_lui)<-c('species','x','y',paste(crop_scen_names[[crop]],lui_scen_names[[lui]],sep=''))
  
  ### loop over species 
    for (s in 1:length(species_names)){
      sp_names<-species_names[[s]]
      sp <- strsplit(sp_names, "_")[[1]][[4]]
      base_dir_out<-paste(path,path_sdm,sp_names,'/',sep='')
  
      # load and crop habitat suitability map
      file<-paste(base_dir_out,'predictions/ensemble/majority_pa/',sp_names,'_ensemble.tif',sep='')
    
      if (file.exists(file)){
      # read eval files an projections
        eval_files <-
          list.files(
            paste0(base_dir_out,"evaluation/"),
            full.names = TRUE,
            recursive = TRUE,
            pattern = "*eval")
      
        evals_out <- lapply(eval_files, get_eval, threshold = "tss")
        eval_df <- do.call(rbind, evals_out)
        eval_df$sp_name <- as.character(eval_df$sp_name)
      
        preds <- list.files(paste0(base_dir_out, "predictions"), full.names = TRUE, recursive = TRUE)
        preds <- preds[!grepl("/ensemble/", preds)]
        preds_r<-lapply(preds, raster)
        preds_c<-lapply(preds_r,crop,extent)
        preds_m<-lapply(preds_c,mask,mask_SA)
      
        # select group-specific LUI maps and apply to projections
        luimaps<- raster::subset(files_LUI,grep(pattern=paste('LUI_map_X2020_ssp2_crop',crop_scen_names[[crop]],sp,sep='_'),names(files_LUI),value=T))

        preds_lui<-list()
        for (i in 1:3){
          preds_lui[[i]] <- preds_m[[i]]*luimaps
        }
      
        ens_preds <- stack(preds_lui) > eval_df$threshold
        ens_pa <- sum(ens_preds)
        ens_pa[ens_pa < round(raster::nlayers(stack(preds_lui)))] <- 0
        ens_pa[ens_pa >= round(raster::nlayers(stack(preds_lui)))] <- 1
        HS_maps_pres_cr<-crop(ens_pa,extent)
        HS_maps_pres_c<-mask(HS_maps_pres_cr, mask_SA)
      
        if (cellStats(HS_maps_pres_c,'sum')>0){
          xy<-as.data.frame(HS_maps_pres_c,xy=TRUE) # raster to data frame
          pres_cells<-as.data.frame(c(sp_names,subset(xy,xy[,3]==1)))
          colnames(pres_cells)<-colnames(df_climate_lui)
          l<-list(df_climate_lui,pres_cells)
          df_climate_lui<-data.table::rbindlist(l,use.names=TRUE)
        }  
        
        else {
          print(paste(sp_names, 'not present in baseline',sep=' '))
        }
      } # if file exists closed
    
      else {
        print(paste(sp_names, ' not available',sep = ''))
      }
    } # species loop closed
  
    outfile_LUI<-paste(path,'11_SDM_output/SA_data_ssp2_',crop_scen_names[[crop]],'.RDS',sep='')# data frame climatic habitat suitability
  
    df_climate_lui<-df_climate_lui[rowSums(is.na(df_climate_lui[,2:4])) < 3,]
    saveRDS(df_climate_lui,outfile_LUI[[crop]])
    print(paste(crop_scen_names[[crop]],lui_scen_names[[lui]], 'complete',sep=' '))
    
  } #crop intensity scenario closed
} # lui (year and ssp) closed

# END

