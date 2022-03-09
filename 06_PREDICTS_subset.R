########################################################################
# script to calculate effect of land use type and intensity on
# animal abundance. This will be used as factor to manipulate SDM output
#
# author: Vivienne Groner
# latest update: 09.03.2022
#########################################################################

webFile <- url("https://timnewbold.github.io/predicts_database.rds?dl=1")
predicts <- readRDS(webFile)
#predicts <- readRDS("predicts_database.rds")
#remotes::install_github("timnewbold/predicts-demo",subdir="predictsFunctions")
#remotes::install_github("timnewbold/StatisticalModels") 
library(predictsFunctions)
library(StatisticalModels)
library("dplyr")                    

predicts <- predictsFunctions::CorrectSamplingEffort(diversity = predicts)
predicts <- predictsFunctions::MergeSites(diversity = predicts,silent = TRUE)

# select unique data for columns:
# 1:20 site information
# 30 "Predominant_land_use" 
# 32 "Use_intensity"
# 37 38 Longitude and latitude
# 58 Phylum 59 Class
# 64 "Best_guess_binomial" 
# 66 "Measurement"  
predicts1 <- unique(predicts[,c(1:20,30,32,37,38,58,59,64,66)])
colnames(predicts1)

# include only data with abundance measurements
predicts2 <- predicts1[(predicts1$Diversity_metric=='abundance'),]

# selection of groups
predicts_arthropoda <- droplevels(predicts2[(predicts2$Phylum=="Arthropoda"),])
predicts_birds <- droplevels(predicts2[(predicts2$Class=="Aves"),])
predicts_mammalia <- droplevels(predicts2[(predicts2$Class=="Mammalia"),])
predicts_reptilia <- droplevels(predicts2[(predicts2$Class=="Reptilia"),])
predicts_amphibia <- droplevels(predicts2[(predicts2$Class=="Amphibia"),])
predicts_gastro <- droplevels(predicts2[(predicts2$Class=="Gastropoda"),])

predicts_animals <- list(predicts_arthropoda,predicts_birds, predicts_mammalia,predicts_amphibia , predicts_reptilia,predicts_gastro)
group <- c('arthropoda', 'birds', 'mammalia', 'amphibia', 'reptilia','gastropoda')
str(predicts_animals)

# start data.frame for factors and records
df <- as.data.frame(matrix(ncol=1))
colnames(df) <- 'LUI_class'
head(df)

for (taxa in 1:length(predicts_animals)){

  # create a site-level data-frame
  sites <- predictsFunctions::SiteMetrics(diversity = predicts_animals[[taxa]],
                                        extra.cols = c("Longitude","Latitude",
                                                       "Predominant_land_use",
                                                       "Use_intensity","SS",
                                                       "SSB","SSBS","Best_guess_binomial"),
                                        srEstimators = NULL)


  #Arranging land-use and intensity classification (from Tim's past code)
  sites$LandUse <- paste(sites$Predominant_land_use)
  sites$LandUse[which(sites$LandUse=="Primary vegetation")] <- "Primary Vegetation"
  sites$LandUse[which(sites$LandUse=="Secondary vegetation (indeterminate age)")] <- NA
  sites$LandUse[which(sites$LandUse=="Cannot decide")] <- NA
  sites$LandUse <- factor(sites$LandUse)
  sites$LandUse <- relevel(sites$LandUse,ref="Primary Vegetation")

  sites$UseIntensity <- paste(sites$Use_intensity)
  sites$UseIntensity[which(sites$Use_intensity=="Cannot decide")] <- NA
  sites$UseIntensity <- factor(sites$UseIntensity)
  sites$UseIntensity <- relevel(sites$UseIntensity,ref="Minimal use")

  #Pasting together land use and land use intensity into a new column
  sites$UI <- paste(sites$LandUse,sites$UseIntensity)
  sites$UI[grep("NA",sites$UI)]<-NA

  sites$UI[which(sites$UI=="Young secondary vegetation Intense use")] <- "Young secondary vegetation Light use"
  sites$UI[which(sites$UI=="Intermediate secondary vegetation Intense use")] <- "Intermediate secondary vegetation Light use"
  sites$UI[which(sites$UI=="Mature secondary vegetation Intense use")] <- "Mature secondary vegetation Light use"
  sites$UI <- factor(sites$UI)
  sites$UI <- relevel(sites$UI,ref="Primary Vegetation Minimal use")
  sites <- sites[!is.na(sites$UI),]

  # group data to calculate mean abundance
  sites <- sites %>%
    group_by(SS,SSB,SSBS,Best_guess_binomial,LandUse,UseIntensity,UI,Longitude, Latitude) %>%
    summarise(Av.abundance=mean(Total_abundance))

  sites$LogAbund <- log((sites$Av.abundance)+1)
  sites$LogAbund[grep("NA",sites$LogAbund)] <- NA
  sites$LogAbund[grep("-Inf",sites$LogAbund)] <- NA
  sites$LogAbund[grep("NaN",sites$LogAbund)] <- NA
  sites <- sites[!is.na(sites$LogAbund),]

  sites$Best_guess_binomial[sites$Best_guess_binomial==""] <- NA
  sites <- sites[!is.na(sites$Best_guess_binomial),]

  # number of data points per land use type
  records <- as.data.frame(cbind(levels(sites$UI),tabulate(sites$UI)))
  colnames(records) <- c('LUI_class', paste(group[[taxa]],'records',sep='_'))

  ### model for land use only
  richMod2 <- StatisticalModels::GLMER(modelData = sites,
                                     responseVar = "LogAbund",
                                     fitFamily = "gaussian",
                                     fixedStruct = "UI",
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     saveVars = c("Longitude","Latitude"),
                                     REM=TRUE)


  richModNull <- StatisticalModels::GLMER(modelData = sites,
                                        responseVar = "LogAbund",
                                        fitFamily = "gaussian",
                                        fixedStruct = "1",
                                        randomStruct = "(1|SS)+(1|SSB)",
                                        saveVars = c("Longitude","Latitude"))



  summary.LM(richMod2)
  summary.LM(richModNull)
  richMod2$model
  richModNull$model
  AIC(richMod2$model,richModNull$model)

  anova(richMod2$model,richModNull$model)

  StatisticalModels::R2GLMER(richMod2$model)

  # plot
  PlotGLMERFactor(model = richMod2$model,data = richMod2$data,
                responseVar = "LogAbund",logLink = "e",
                catEffects = "UI",xtext.srt = 35,
                #order = c(1,8),
                main=group[[taxa]],
                params = list(mar = c(3, 3.3, 1, 1)))

  # write data frame UI and model output
  LUI_class_taxa <- data.frame(gsub('UI','',dimnames(richMod2$model@pp$X)[[2]]))
  LUI_class_taxa$data<-richMod2$model@beta
  colnames(LUI_class_taxa) <- c('LUI_class',group[[taxa]])

  # merge records and LUI and model output
  df_list <- list(LUI_class_taxa,records,df)
  df <- Reduce(function(x, y) merge(x, y,by=c('LUI_class'),all=TRUE), df_list)

}

df
write.csv(df,file=paste(my_path,'02_environmental_data/PREDICTS/LUI_groups.csv',sep=''))

# END
