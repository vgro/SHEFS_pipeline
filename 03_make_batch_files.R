################################################################
# write batch file to submit multiple parallel jobs on cluster 
# xxxxxx is your ID
#
# author: Vivienne Groner
# latest update: 09.03.2022
################################################################

# read species names
wd <- paste(my_path)
files <- list.files(paste(wd,'04_occurrence_records/GBIF_ready/points/rarefied/',sep=''))
files


batch_test <-
'#!/bin/bash -l
#$ -l h_rt=01:00:0
#$ -l mem=12G
#$ -l tmpfs=4G
#$ -N test_batch
#$ -wd /home/xxxxxx/Scratch/SHEFS
#echo $JOB_ID > /home/xxxxxx/Scratch/SHEFS/x_jobname_$JOB_ID.txt
'

st_f <- 1
end_f <- 25
for (f in 1:72){
  
  group1 <- files[st_f:end_f]
  entry <- list()
  for (i in 1:length(group1)){
    
    entry_name <- gsub(".csv", ".sh", group1[[i]])
    entry[[i]] <- paste('qsub /home/xxxxxx/Scratch/SHEFS/sh_scripts/',entry_name,sep='')
  }
  
  dt <- unlist(entry)
  job <- gsub("test_batch", paste('batch_',st_f,'-',end_f,sep=''), batch_test)
  write(c(batch_test,dt),file=paste(wd,'move_to_myriad/batch_scripts/batch_',st_f,'-',end_f,'.sh',sep=''))
  
  st_f <- st_f+25
  end_f <- end_f+25
  
}

# END

