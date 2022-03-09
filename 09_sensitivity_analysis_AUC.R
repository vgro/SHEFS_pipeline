####################################
## AUC distribution
## author: Vivienne Groner
#
# Last update: 09.03.2022
###################################

# load libraries
library(stringr)
library(ggplot2)


#set environment
path <- my_path
source(paste0(my_path,"05_code_vivienne/pipeline_and_analysis/SA_functions_CBER_analysis.R"))
memory.limit(size=500000000)


# read species names
sp_names <- gsub(".csv", "", list.files(paste0(path,"/04_occurrence_records/GBIF_ready/points/rarefied/"),pattern='csv'))


# load evaluation files and get AUCs
eval_files<-list()

for (s in 1:length(sp_names)){
  
  eval_files[[s]] <-
  list.files(
    paste0(path,"11_SDM_output/data/",sp_names[[s]],'/evaluation/bioclim/'),
    full.names = TRUE,
    recursive = TRUE,
    pattern = "*eval.RDA")
  
}

evals_out <- lapply(unlist(eval_files), get_eval, threshold = "tss")
eval_df <- do.call(rbind, evals_out)


# plot histogram of all AUCs
p<-ggplot(eval_df, aes(x=auc)) + 
  geom_histogram(binwidth=0.05,col='black',fill='grey')
p + scale_color_grey() + theme_bw() 

# END