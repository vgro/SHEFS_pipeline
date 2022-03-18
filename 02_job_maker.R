###################################################### 
# create R scripts and sh scripts for each species 
# to run SDM pipeline on cluster 
#
# author: Vivienne Groner
# date: 18.03.2022
########################################################

# set environment and load dummy files
wd <- paste(my_path)
files<-list.files(paste0(wd,'04_occurrence_records/GBIF_ready/points/rarefied/'))
out_dir_R <- paste0(wd,'05_code/r_jobs/')
out_dir_sh <- paste0(wd,'05_code/sh_scripts/')

dummy_sh <- paste(wd,'05_code/SA_dummy.sh',sep='')
dummy_R <- paste(wd,'05_code/SA_dummy.R',sep='')   


# create R scripts
for(i in 1:length(files)){
  
  ptc <- readLines(dummy_R)  #Read in the template file
  job <- gsub(".csv", "", files[[i]])
  ptc[17] <- gsub("dummy", job, ptc[17])
  file_out <- gsub(".csv", ".R", files[[i]])
  writeLines(ptc, paste(out_dir_R,file_out,sep=''))  #Use this as the file name and write it to the output directory

  }


# bash_scripts to submit all the jobs
for(j in 1:length(files)){#
  
  ptc <- readLines(dummy_sh)  #Read in the template file
  job <- gsub(".csv", "", files[[j]])
  ptc[5] <- gsub("dummy", job, ptc[5])
  ptc[22] <- gsub("dummy", job, ptc[22])
  ptc[23] <- gsub("dummy", job, ptc[23])
  ptc[24] <- gsub("dummy", job, ptc[24])
  ptc[25] <- gsub("dummy", job, ptc[25])
  ptc[26] <- gsub("dummy", job, ptc[26])
  ptc[27] <- gsub("dummy", job, ptc[27])
  ptc[28] <- gsub("dummy", job, ptc[28])
  
  writeLines(ptc, paste0(out_dir_sh,job, ".sh", sep = ""))  #Use this as the file name and write it to the output directory

  }

# END