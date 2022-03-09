##########################################################
## test buffer and density in pseudoabsences selection
## author: Vivienne Groner
#
# Last update: 09.03.2022
##########################################################

# load libraries
library(dplyr)
library(raster)
library(lwgeom)
library(sp)
library(readr)
library(sf)
library(dismo)
library(rgdal)  
library(stringr)
library(ggplot2)


#set environment
path <- my_path
source(paste0(my_path,"05_code_vivienne/pipeline_and_analysis/SA_functions_CBER_analysis.R"))
memory.limit(size=500000000)

if(!dir.exists(paste(path,"/04_occurrence_records/buffer_test/",sep=''))){
  dir.create(paste(path,"/04_occurrence_records/buffer_test/",sep=''))
}

# study area
xmin <- 15
xmax <- 35
ymin <- -35
ymax <- 0

ext_occ <- raster::extent(xmin, xmax, ymin, ymax)
  

# load ecoregions  
ecoreg <-
  sf::st_read(paste(path,"/02_environmental_data/WWF_Ecoregions/wwf_terr_ecos.shp",sep='')) %>% 
  sf::st_crop(.,ext_occ) %>%  ##cropping to the area of interest
  dplyr::select(OBJECTID, ECO_NAME) ##just selecting out the columns we're interested in


# load climate data
climate_path <- paste0(path,'/02_environmental_data/CHELSA')
bio_layers <- list.files(climate_path, pattern = 'tif',full.names=TRUE)
bio_layers_r <-lapply(bio_layers,raster::raster)
bio_layers_s <- raster::stack(bio_layers_r)
env_crop <- raster::crop(bio_layers_s, ext_occ)

bioclim_layers <- c(1, 5, 6, 13, 14)

# get species names
sp_names_in <- gsub(".csv", "", list.files(paste0(path,"/04_occurrence_records/buffer_test/rarefied/"),pattern='csv'))
sp_names <- sample(x=sp_names_in,200)


# sample pseudoabsences over a range of buffer and density values
density_range <- c(100,250,500)
buffer_range <- c(10,25,50)

for (b in 1:length(buffer_range)){ 
  
  for (d in 1:length(density_range)){ 
    
    lapply(X = sp_names, 
       background_sampler_sens, 
       in_dir = paste(path,"/04_occurrence_records/buffer_test/rarefied/",sep=''), 
       out_dir = paste(path,"/04_occurrence_records/buffer_test/pseudo/",sep=''),
       dens_abs = "density", 
       density = density_range[[d]], 
       type = "pseudoabsence", 
       buffer = buffer_range[[b]], 
       polygon = ecoreg)
 
  }

}

# load available pseudoabsences
sp_names <- unique(str_sub(list.files(paste0(path,"/04_occurrence_records/buffer_test/pseudo/"),pattern='csv'),1,-18))

# check  for complete data
for (s in 1:length(sp_names)){
  
  sp_in<-list.files(paste0(path,"/04_occurrence_records/buffer_test/pseudo/"),pattern=sp_names[[s]],full.names=TRUE)
  
  if (length(sp_in)<9){
    
    print(paste(sp_names[[s]],'not complete'),sep=' ')
    
  }
  
}


# run bioclim SDM with different background points
for (b in 1:length(buffer_range)){
  
  for (d in 1:length (density_range)){
    
    lapply(
      X = sp_names,
      fitBC_senstest,
      pres_dir = paste0(path,"/04_occurrence_records/GBIF_ready/environmental/presence/"),
      backg_dir = paste0(path,"/04_occurrence_records/buffer_test/pseudo/"),
      predictor_names = bioclim_layers,
      predictors = env_crop,
      pred_out_dir = paste0(path,"/04_occurrence_records/buffer_test/predictions/"),
      eval_out_dir = paste0(path,"/04_occurrence_records/buffer_test/predictions/"),
      overwrite = TRUE,
      threads = 4,
      eval = TRUE,
      buf=buffer_range[[b]],
      den=density_range[[d]]
    )
  
  }
  
}

#############################################################
# analysis
#############################################################

# check for complete data (projections)
sp_names <- unique(str_sub(list.files(paste0(path,"/04_occurrence_records/buffer_test/predictions/"),pattern='tif'),1,-20))

for (s in 1:length(sp_names)){
  
  sp_in<-list.files(paste0(path,"/04_occurrence_records/buffer_test/predictions/"),pattern=sp_names[[s]],full.names=TRUE)
  
  if (length(sp_in) < 27){
    
    print(paste(sp_names[[s]],'not complete'),sep=' ')
    
  }
  
}


# load evaluation files and get AUCs
eval_files <-
  list.files(
    paste0(path,"/04_occurrence_records/buffer_test/predictions/"),
    full.names = TRUE,
    recursive = TRUE,
    pattern = "*eval.RDS")

evals_out <- lapply(eval_files, get_eval, threshold = "tss")
eval_df <- do.call(rbind, evals_out)


# create repeating test names for buffer and density values for data frame
test_names <- as.vector(replicate(length(sp_names),c('den100_buf10','den250_buf10','den500_buf10',
             'den100_buf25','den250_buf25','den500_buf25',
             'den100_buf50','den250_buf50','den500_buf50')))


# combine evaluation output and test names             
eval_df_test <- cbind(test_names,eval_df)
head(eval_df_test)


# grep scenario save histogram of AUCs
for (i in 1:3){
  
  for (j in 1:3){
    
    eval_df_single<-subset(eval_df_test,eval_df_test$test_names==paste0('den',density_range[[i]],'_buf',buffer_range[[j]]))
    eval_df_single
    
    ggplot(eval_df_single, aes(x=auc)) + 
      geom_histogram(binwidth=0.05,col='black',fill='grey') + 
      ylim(0, 50)+
      scale_color_grey() + theme_bw() +
      annotate(geom="text", x=0.5, y=49, label=paste0('density = ',density_range[[i]]))+
      annotate(geom="text", x=0.5, y=46, label=paste0('buffer  =    ',buffer_range[[j]]))
    
    ggsave(paste0(path,"/04_occurrence_records/buffer_test/AUC_hist_den",density_range[[i]],'_buf',buffer_range[[j]],".png"))
  }
  
} 
  

# calculate mean and standard deviation of auc over all species
eval_grouping<-list()
for (i in 1:9){
  
  eval_grouping[[i]] <- data.frame(test_names[[i]],
                                 median(subset(eval_df_test$auc,eval_df_test$test_names==test_names[[i]])),
                                 mean(subset(eval_df_test$auc,eval_df_test$test_names==test_names[[i]])),
                                 sd(subset(eval_df_test$auc,eval_df_test$test_names==test_names[[i]])))
  
  colnames( eval_grouping[[i]]) <- c('test','median','mean','sd')
  
}

eval_means <- do.call(rbind,eval_grouping)
head(eval_means)


# plot mean and sd over all species
ggplot(eval_means, aes(x=test, y=median)) +
  geom_point(position=position_dodge(width = .2), stat="identity",
           colour='black') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2)+
  theme_bw()

ggsave(paste0(path,"/04_occurrence_records/buffer_test/AUC_all_species_all_denbuf.png"))


# compare k-fold aucs for each species individually
evals_out_auc <- lapply(eval_files, get_eval_auc)
eval_df_auc <- do.call(rbind, evals_out_auc)
eval_df_auc <- cbind(test_names,eval_df_auc)
head(eval_df_auc)


# histogram of all mean aucs
p<-ggplot(eval_df_auc, aes(x=aucs)) + 
  geom_histogram(binwidth=0.05,col='black',fill='grey')
p + scale_color_grey() + theme_bw() 


# effect of buffer and density vs k-fold
df_both_auc <- data.frame(matrix(ncol=3))
colnames(df_both_auc) <- c('kfold','buffer', 'density')

for (s in 1:length(sp_names)){
  
  data1 <- subset(eval_df_auc,eval_df_auc$sp_name == sp_names[[s]])
  df_test <- data.frame(matrix(ncol=3))
  colnames(df_test) <- c('kfold','buffer','density')
  
  df_test$kfold <- mean(data1$sd)
  df_test$buffer <- sd(c(data1$aucs[1:3],data1$aucs[4:6],data1$aucs[7:9]))
  df_test$density <- sd(c(data1$aucs[c(1,4,7)],data1$aucs[c(2,5,8)],data1$aucs[c(3,6,9)]))
  
  l <- list(df_both_auc,df_test)
  df_both_auc <- data.table::rbindlist(l,use.names=TRUE)
  
}
 
df_both_auc <- df_both_auc[rowSums(is.na(df_both_auc[,1:2])) < 2,]
head(df_both_auc)

boxplot(df_both_auc)

# END
    
  