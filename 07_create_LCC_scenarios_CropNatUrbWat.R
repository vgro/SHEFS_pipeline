#######################################################
# script to create land use and intensity suitability
# and future land cover scenarios
# natural=1-urban-cropland-water
# 
# author: Vivienne Groner
# latest update: 09.03.2022
#######################################################

# set environment and load libraries
library(raster)
library(sp)
library(rgdal)
library(RColorBrewer)

path_in  <- paste0(my_path,'02_environmental_data/')
path_out <- paste0(my_path,'02_environmental_data/SSP_landcover/')


# reference for resampling
ref_file <- raster(paste0(my_path,'11_SDM_output/data/Acherontia_atropos_GBIF_arthropoda/predictions/ensemble/majority_pa/Acherontia_atropos_GBIF_arthropoda_ensemble.tif'))

extent <- c(15, 34, -35, -21)
years <- c(2020,2080)
ssp <- c('ssp2','ssp5')
rcp <- c('RCP45','RCP85')


# mask SA shape file
mask_SA_in<-readOGR(paste(path_in,'Igismap/SouthAfrica_Bound.shp',sep=''))
mask_SA<-crop(mask_SA_in,extent)
plot(mask_SA)

# load cropland maps
# Cao, B.et al.: A 1 km global cropland dataset from 10 000 BCE to 2100 CE, 
# Earth Syst. Sci. Data, 13, 5403–5421, https://doi.org/10.5194/essd-13-5403-2021, 2021. 

cropland <- list()
i=1
for (y in 1:length(years)){
  
  for (s in 1:length(ssp)){
    
    cropland_in <- raster(paste(path_in,'Cao_2021_cropland/globalCropland_',years[[y]],'CE_',ssp[[s]],'_',rcp[[s]],'.tif',sep=''))
    cropland[[i]] <- crop(cropland_in,extent)
    names(cropland[[i]]) <- paste('cropland',years[[y]],ssp[[s]],sep='_')
    i=i+1
    
  }
  
}

# load urban area 
# Gao, J. and M. Pesaresi. 2021. Global 1-km Downscaled Urban Land Extent 
# Projection and Base Year Grids by SSP Scenarios, 2000-2100. Palisades, NY: 
# NASA Socioeconomic Data and Applications Center (SEDAC). https://doi.org/10.7927/1z4r-ez63. 

urban<-list()
i=1
for (y in 1:length(years)){
  
  for (s in 1:length(ssp)){
    
    urban_in <- raster(paste(path_in,'Gao_2021_urban/UrbanFraction_1km_GEOTIFF_Projections_SSPs1-5_2010-2100_v1/',ssp[[s]],'_',years[[y]],'.tif',sep=''))
    urban[[i]] <- crop(urban_in,extent)
    names(urban[[i]]) <- paste('urban',years[[y]],ssp[[s]],sep='_')
    i=i+1
    
  }
  
}


# load water area
# South African National Land Cover (SANLC) 2020
water_in <- raster(paste(path_in,'HighRes_SA_landcover_map/SANLC_2020_GEOTIFF/SA_NLC_resampled/Resamp_water.tif',sep=''))
water <- crop(water_in,extent)


# (semi)natural habitat as 1-cropland-urban-water
natural_matrix1 <- raster(matrix(1, nrow=1536,ncol=2004))
extent(natural_matrix1) <- extent
crs(natural_matrix1) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
natural_matrix <- resample(natural_matrix1,cropland[[1]])
natural_matrix

natural <- list()
i=1
for (y in 1:length(years)){
  
  for (s in 1:length(ssp)){
    
    natural[[i]] <- natural_matrix-cropland[[i]]-urban[[i]]-water
    natural[[i]][natural[[i]]<0] <- 0
    names(natural[[i]]) <- paste('natural',years[[y]],ssp[[s]],sep='_')
    i=i+1
    
  }
  
}


# plot for eyeballing
par(mfrow=c(2,2))
plot(mask(natural[[1]],mask_SA),main='natural')
plot(mask(cropland[[1]],mask_SA), main='cropland')
plot(mask(urban[[1]],mask_SA), main='urban')
plot(mask(water,mask_SA),main='water')
dev.off()

# habitat suitability factors
# urban group dependent
# cropland group dependent
# natural 1
# water 0
# baseline is ssp245  !!!

group <- c('arthropoda', 'birds', 'mammalia', 'amphibia', 'reptilia','gastropoda')
crop_intense_f <- c(0.1,1,1,0.5,0.5,0.5)
urban_f <- c(0.4,1,1,0,0,0)


# plot settings
palette = brewer.pal(11,'YlGn')
data_frame <- data.frame(name = group, # write group in plot
                         latitude = rep(-22,6),
                         longitude = rep(19.2,6))


for( g in 1:length(group)){
  
  for (y in 1:length(natural)){
    
    LUI_map_intense <- crop_intense_f[[g]]*cropland[[y]]+urban_f[[g]]*urban[[y]]+natural[[y]]
    LUI_map_intense[LUI_map_intense>1] <- 1
    LUI_map_intense <- resample(LUI_map_intense,ref_file)
    names(LUI_map_intense) <- paste(substr(names(natural[[y]]),9,nchar(names(natural[[y]]))),'_crop_intense_',group[[g]],sep='')
    writeRaster(LUI_map_intense,paste(path_out,'LUI_map_',names(LUI_map_intense),'.tif',sep=''),overwrite=TRUE)
    print(names(LUI_map_intense))
    
    # this part saves the plot
    data_frame1 <- data.frame(name = gsub('_',' ',substr(names(natural[[y]]),9,nchar(names(natural[[y]])))), # write scenario in plot
                              latitude = -21,
                              longitude = 19)
    
    
    tiff(paste(path_out,'plots/LUI_map_',names(LUI_map_intense),'.png',sep='')) 
    plot(mask(LUI_map_intense,mask_SA),xlab="longitude", ylab="latitude",col=palette,
         xlim=c(16.4,34),ylim=c(-35,-21.5))
    text(x = data_frame$longitude, y = data_frame$latitude, 
         data_frame$name[[g]], pos = 1, col = "black")
    text(x = data_frame1$longitude, y = data_frame1$latitude, 
         data_frame1$name, pos = 1, col = "black")
    dev.off()

  }
  
}

# END
