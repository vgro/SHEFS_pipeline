#!/bin/bash -l
#$ -l h_rt=08:00:0
#$ -l mem=12G
#$ -l tmpfs=4G
#$ -N dummy
#$ -wd /home/xxxxxx/Scratch/SHEFS


#echo $JOB_ID > /home/xxxxxx/Scratch/SHEFS/x_jobname_$JOB_ID.txt

cp -a /home/xxxxxx/Scratch/SHEFS/R_Library/ $TMPDIR
cp -a /home/xxxxxx/Scratch/SHEFS/CHELSA/ $TMPDIR
cp -a /home/xxxxxx/Scratch/SHEFS/CHELSA_future/ $TMPDIR
cp /home/xxxxxx/Scratch/SHEFS/SA_functions_CBER.R $TMPDIR

cd $TMPDIR

module unload compilers
module unload mpi
module load r/recommended

R --no-save < /home/xxxxxx/Scratch/SHEFS/r_jobs/dummy.R > test.out
tar zcvf ~/Scratch/SHEFS/dummy.tar.gz ~/Scratch/SHEFS/dummy/
rm -r ~/Scratch/SHEFS/dummy/
rm ~/Scratch/SHEFS/dummy.e$JOB_ID
rm ~/Scratch/SHEFS/dummy.o$JOB_ID
rm ~/Scratch/SHEFS/dummy_log.txt



