###############################################################################
# script to create land cover suitability map for SHEFS joint project
#
# 1. resample SA_NLC_2020 to 1km resolution and calculate % of LCT per grid cell
#   (reclassified in ArcMap)
# 2. calculate suitability based on PREDICTS database (effects on abundance, 
# following Newbold et al. 2015)
#
# author: Vivienne Groner
# latest update: 09.03.2022
##############################################################################

library(raster)
path <- paste0(my_path,'02_environmental_data/')

# load raster land cover types
infile <- raster(paste(path,'HighRes_SA_landcover_map/SANLC_2020_GEOTIFF/SA_NLC_2020_PREDICTS.tif',sep=''))
infile

# climate raster as reference for resampling
reference_file <- raster(paste(path,'CHELSA/CHELSA_bio10_01.tif',sep=''))
reference_file

ref_c <- crop(reference_file,infile)
ref_c

# aggregate to 1 km and calculate percentage of land cover class per grid cell
cov_pct <- lapply(unique(infile), function(land_class) {
  aggregate(infile, fact=46, fun=function(vals, na.rm) {
    sum(vals==land_class, na.rm=na.rm)/length(vals)
  })
})

# only the PREDICTS land use and intensity types that are represented in the data set
names(cov_pct) <- c('bare','primary-minimal','primary-light','secondary-minimal','plantation-minimal','plantation-light',
                  'plantation-intense','cropland-light','cropland-intense','urban-minimal','urban-light',
                  'urban-intense','water','none')


# resample to reference (climate) raster
resampled <- list()
for (i in 1:13){
  
  resampled[[i]]<-resample(cov_pct[[i]],ref_c, method='bilinear')
  print(names(cov_pct)[[i]])
  
}

# check if dimensions are correct
ref_c
resampled[[1]]

# save raster for each land cover type
for (i in 1:13){
  
  writeRaster(resampled[[i]],file=paste(path, 'HighRes_SA_landcover_map/SANLC_2020_GEOTIFF/PREDICTS_resampled/Resamp_',names(cov_pct)[[i]],'PREDICTS.tif',sep=''))

}

# END
