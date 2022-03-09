.libPaths(c('/lustre/scratch/scratch/xxxxxx/SHEFS/R_Library',.libPaths()))

library(dplyr)
library(dismo)
library(here)
library(lwgeom)
library(raster)
library(readr)
library(sf)
library(sp)
library(stringr)
library(vroom)
library(ggplot2)

source("/lustre/scratch/scratch/xxxxxx/SHEFS/SA_functions_CBER.R")

sp_names <- 'dummy'
base_dir_in <- "/lustre/scratch/scratch/xxxxxx/SHEFS/"

if(!dir.exists(paste0(base_dir_in,sp_names,"/"))){
  dir.create(paste0(base_dir_in,sp_names,"/"), recursive = TRUE)
}
base_dir_out<-paste(base_dir_in,sp_names,'/',sep='')


###Data loading and prepping
climate_path <- paste(base_dir_in,'CHELSA/',sep='')
climate_path_future <- c(paste0(base_dir_in,'CHELSA_future/2061-2080/rcp45/ACCESS1-0/'),
                       paste0(base_dir_in,'CHELSA_future/2061-2080/rcp45/CNRM-CM5/'),
                       paste0(base_dir_in,'CHELSA_future/2061-2080/rcp45/GFDL-ESM2G/'),
                       paste0(base_dir_in,'CHELSA_future/2061-2080/rcp45/HadGEM2-CC/'),
                       paste0(base_dir_in,'CHELSA_future/2061-2080/rcp45/MPI-ESM-MR/'),
                       paste0(base_dir_in,'CHELSA_future/2061-2080/rcp85/ACCESS1-0/'),
                       paste0(base_dir_in,'CHELSA_future/2061-2080/rcp85/CNRM-CM5/'),
                       paste0(base_dir_in,'CHELSA_future/2061-2080/rcp85/GFDL-ESM2G/'),
                       paste0(base_dir_in,'CHELSA_future/2061-2080/rcp85/HadGEM2-CC/'),
                       paste0(base_dir_in,'CHELSA_future/2061-2080/rcp85/MPI-ESM-MR/')
                       )

climate_scen_name <- c('rcp45_ACCESS1-0','rcp45_CNRM-CM5','rcp45_GFDL-ESM2G','rcp45_HadGEM2-CC','rcp45_MPI-ESM-MR',
                     'rcp85_ACCESS1-0','rcp85_CNRM-CM5','rcp85_GFDL-ESM2G','rcp85_HadGEM2-CC','rcp85_MPI-ESM-MR')
xmin <- 0
xmax <- 45
ymin <- -45
ymax <- 0

min_occ <- 15
bioclim_layers <- c(1, 5, 6, 13, 14)#the number of the bioclim layers to be included as environmental variables - https://worldclim.org/bioclim
ext_occ <- extent(xmin, xmax, ymin, ymax)


# predictors for calibration present day
bio_layers <- list.files(climate_path, pattern = 'tif')
bio_layers1<-paste(climate_path, bio_layers, sep="") 
bio_layers2 <-lapply(bio_layers1,raster)
env_layers <- raster::stack(bio_layers2)
env_crop <- raster::crop(env_layers, ext_occ)


# predictors for future climate
env_crop_future <- list()
for (fut in 1:length(climate_scen_name)){
  
  bio_layers_future <- list.files(climate_path_future[[fut]], pattern = 'tif')
  bio_layers1_future <- paste0(climate_path_future[[fut]], bio_layers_future) 
  bio_layers2_future <- lapply(bio_layers1_future,raster)
  env_layers_future <- raster::stack(bio_layers2_future)
  env_crop_future[[fut]] <- raster::crop(env_layers_future, ext_occ)
  names(env_crop_future[[fut]]) <- c('CHELSA_bio10_01','CHELSA_bio10_13','CHELSA_bio10_14','CHELSA_bio10_05','CHELSA_bio10_06')

}


# Create folders for predictions
if(!dir.exists(paste0(base_dir_out,"predictions/bioclim/"))){
  dir.create(paste0(base_dir_out,"predictions/bioclim"), recursive = TRUE)
}

if(!dir.exists(paste0(base_dir_out, "predictions/glm/"))){
  dir.create(paste0(base_dir_out,"predictions/glm"), recursive = TRUE)
}

if(!dir.exists(paste0(base_dir_out,"predictions/rf/"))){
  dir.create(paste0(base_dir_out,"predictions/rf"), recursive = TRUE)
}
#
if(!dir.exists(paste0(base_dir_out,"evaluation/bioclim/"))){
  dir.create(paste0(base_dir_out,"evaluation/bioclim"), recursive = TRUE)
}

if(!dir.exists(paste0(base_dir_out, "evaluation/glm/"))){
  dir.create(paste0(base_dir_out,"evaluation/glm"), recursive = TRUE)
}

if(!dir.exists(paste0(base_dir_out,"evaluation/rf/"))){
  dir.create(paste0(base_dir_out,"evaluation/rf"), recursive = TRUE)
}


for (fut in 1:length(climate_scen_name)){
  
  if(!dir.exists(paste0(base_dir_out,"future/",climate_scen_name[[fut]],"/predictions/bioclim/"))){
    dir.create(paste0(base_dir_out,"future/",climate_scen_name[[fut]],"/predictions/bioclim/"), recursive = TRUE)
  }
  
  if(!dir.exists(paste0(base_dir_out, "future/",climate_scen_name[[fut]],"/predictions/glm/"))){
    dir.create(paste0(base_dir_out,"future/",climate_scen_name[[fut]],"/predictions/glm/"), recursive = TRUE)
  }
  
  if(!dir.exists(paste0(base_dir_out,"future/",climate_scen_name[[fut]],"/predictions/rf/"))){
    dir.create(paste0(base_dir_out,"future/",climate_scen_name[[fut]],"/predictions/rf/"), recursive = TRUE)
  }
  
}


#### Fit Bioclim Models
# present
fitBC(
  X = sp_names,
  pres_dir = paste0(base_dir_in,"environmental/presence/"),
  backg_dir = paste0(base_dir_in, "environmental/pseudoabsence/"),
  predictor_names = bioclim_layers,
  predictors = env_crop,
  pred_out_dir = paste0(base_dir_out,"predictions/bioclim/"),
  eval_out_dir = paste0(base_dir_out,"evaluation/bioclim/"),
  overwrite = TRUE,
  threads = 4,
  eval = TRUE
)


# future
for (fut in 1:length(climate_scen_name)){
  
  fitBC_future(
    X = sp_names,
    pres_dir = paste0(base_dir_in,"environmental/presence/"),
    backg_dir = paste0(base_dir_in, "environmental/pseudoabsence/"),
    predictor_names = bioclim_layers,
    predictors_future = env_crop_future[[fut]],
    pred_out_dir_future = paste0(base_dir_out,"future/",climate_scen_name[[fut]],"/predictions/bioclim/"),
    eval_out_dir = paste0(base_dir_out,"evaluation/bioclim/"),
    overwrite = TRUE,
    threads = 4,
    eval = TRUE
  )
  
}


#### Fit GLM Models
# present
fitGLM(
  X = sp_names,
  pres_dir = paste0(base_dir_in,"environmental/presence/"),
  backg_dir = paste0(base_dir_in, "environmental/pseudoabsence/"),
  predictor_names = bioclim_layers,  
  predictors = env_crop,
  pred_out_dir =  paste0(base_dir_out,"predictions/glm/"),
  eval_out_dir =  paste0(base_dir_out,"evaluation/glm/"),
  overwrite = TRUE,
  threads = 4,
  eval = TRUE
)


# future
for (fut in 1:length(climate_scen_name)){
  
  fitGLM_future(
    X = sp_names,
    pres_dir = paste0(base_dir_in,"environmental/presence/"),
    backg_dir = paste0(base_dir_in, "environmental/pseudoabsence/"),
    predictor_names = bioclim_layers,  
    predictors_future = env_crop_future[[fut]],
    pred_out_dir_future = paste0(base_dir_out,"future/",climate_scen_name[[fut]],"/predictions/glm/"),
    eval_out_dir =  paste0(base_dir_out,"evaluation/glm/"),
    overwrite = TRUE,
    threads = 4,
    eval = TRUE
  )
  
}


#### Fit random forest models
# present
fitRF(
 X = sp_names,
 pres_dir =  paste0(base_dir_in,"environmental/presence/"),
 backg_dir =  paste0(base_dir_in,"environmental/pseudoabsence/"),
 predictor_names = bioclim_layers,
 predictors = env_crop,
 pred_out_dir =  paste0(base_dir_out,"predictions/rf/"),
 eval_out_dir =  paste0(base_dir_out,"evaluation/rf/"),
 overwrite = TRUE,
 threads = 4,
 eval = TRUE
)


# future
for (fut in 1:length(climate_scen_name)){
  
  fitRF_future(
    X = sp_names,
    pres_dir =  paste0(base_dir_in,"environmental/presence/"),
    backg_dir =  paste0(base_dir_in,"environmental/pseudoabsence/"),
    predictor_names = bioclim_layers,
    predictors_future = env_crop_future[[fut]],
    pred_out_dir_future = paste0(base_dir_out,"future/",climate_scen_name[[fut]],"/predictions/rf/"),
    eval_out_dir =  paste0(base_dir_out,"evaluation/rf/"),
    overwrite = TRUE,
    threads = 4,
    eval = TRUE
  )  
  
}


#### Get evaluation and AUCs out
eval_files <-
  list.files(
    paste0(base_dir_out,"evaluation/"),
    full.names = TRUE,
    recursive = TRUE,
    pattern = "*eval")

evals_out <- lapply(eval_files, get_eval, threshold = "tss")
eval_df <- do.call(rbind, evals_out)
eval_df$sp_name <- as.character(eval_df$sp_name)

####
if(!dir.exists( paste0(base_dir_out,"predictions/ensemble/majority_pa"))){
  dir.create(paste0(base_dir_out,"predictions/ensemble/majority_pa"), recursive = TRUE)
}

if(!dir.exists(paste0(base_dir_out,"predictions/ensemble/weighted"))){
  dir.create(paste0(base_dir_out,"predictions/ensemble/weighted"), recursive = TRUE)
}


####
for (fut in 1:length(climate_scen_name)){
  
  if(!dir.exists( paste0(base_dir_out,"future/",climate_scen_name[[fut]],"/predictions/ensemble/majority_pa"))){
  dir.create(paste0(base_dir_out,"future/",climate_scen_name[[fut]],"/predictions/ensemble/majority_pa"), recursive = TRUE)
  }

  if(!dir.exists(paste0(base_dir_out,"future/",climate_scen_name[[fut]],"/predictions/ensemble/weighted"))){
  dir.create(paste0(base_dir_out,"future/",climate_scen_name[[fut]],"/predictions/ensemble/weighted"), recursive = TRUE)
  }
  
}

preds <- list.files(paste0(base_dir_out, "predictions"), full.names = TRUE, recursive = TRUE)
preds <- preds[!grepl("/ensemble/", preds)]

lapply(X = sp_names, 
       FUN = ensemble_model, 
       eval_df = eval_df, 
       preds = preds, 
       out_dir = paste0(base_dir_out, "predictions/ensemble/"), 
       method = "majority_pa")

lapply(X = sp_names, 
       FUN = ensemble_model, 
       eval_df = eval_df, 
       preds = preds, 
       out_dir = paste0(base_dir_out, "predictions/ensemble/"), 
       method = "weighted")


# future 
for (fut in 1:length(climate_scen_name)){
  
  predsf <- list.files(paste0(base_dir_out, "future/",climate_scen_name[[fut]],"/predictions"), full.names = TRUE, recursive = TRUE)
  predsf <- predsf[!grepl("/ensemble/", predsf)]

 lapply(X = sp_names, 
        FUN = ensemble_model_f, 
        eval_df = eval_df, 
        preds = predsf, 
        out_dir = paste0(base_dir_out, "future/",climate_scen_name[[fut]],"/predictions/ensemble/"), 
        method = "majority_pa")
 
 lapply(X = sp_names, 
        FUN = ensemble_model_f, 
        eval_df = eval_df, 
        preds = predsf, 
        out_dir = paste0(base_dir_out, "future/",climate_scen_name[[fut]],"/predictions/ensemble/"), 
        method = "weighted")
 
}


########Plots for eyeballing
if(!dir.exists(paste0(base_dir_out, "plots/ensemble/majority_pa"))){
  dir.create(paste0(base_dir_out, "plots/ensemble/majority_pa"), recursive = TRUE)
}

for ( fut in 1:length(climate_scen_name)){
  
  if(!dir.exists(paste0(base_dir_out, "plots/future/",climate_scen_name[[fut]],"/ensemble/majority_pa"))){
  dir.create(paste0(base_dir_out, "plots/future/",climate_scen_name[[fut]],"/ensemble/majority_pa"), recursive = TRUE)
  }
  
}

ggplot_out(X = sp_names, 
           points_dir = paste0(base_dir_in, "points/rarefied/"), 
           rast_dir = paste0(base_dir_out,"predictions/ensemble/majority_pa/"), 
           out_dir = paste0(base_dir_out, "plots/ensemble/majority_pa/"))


for (fut in 1:length(climate_scen_name)){
  
  lapply(X = sp_names, 
         FUN = ggplot_out, 
         points_dir = paste0(base_dir_in, "points/rarefied/"), 
         rast_dir = paste0(base_dir_out,"future/",climate_scen_name[[fut]],"/predictions/ensemble/majority_pa/"), 
         out_dir = paste0(base_dir_out, "plots/future/",climate_scen_name[[fut]],"/ensemble/majority_pa/"))
  
}


##### end #######
