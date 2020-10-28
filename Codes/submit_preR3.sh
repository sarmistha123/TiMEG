#!/bin/bash
#----------------------------------------------------------
# Job name
#PBS -N testjob 

# Queue name
#PBS -q parallel

# Name of stdout output file
##PBS -o job.out

# Total number of nodes and MPI tasks/node requested
#PBS -l select=1:ncpus=1
##PBS -l select=4:mpiprocs=40



# Run time (hh:mm:ss) - 1.5 hours
# PBS -l walltime=01:30:00
#----------------------------------------------------------

# Change to submission directory
cd $PBS_O_WORKDIR

#export OMP_NUM_THREADS=40
# Launch MPI-based executable
#prun R CMD BATCH first.R
#bbb=`ls -trl | ls  fileid* | tail -1 | awk '{print $NF}' | cat`
#bbb=`ls | head | grep fileid.*`
bbb=`ls -trl | ls  fileid* | head -1 | awk '{print $NF}' | cat`
name="$(cut -d. -f3-6 <<< `echo $bbb`)"
mv $bbb $name
#echo $name
Rscript TiMEG.R $name
#prun hw_omp_mpbrid.exe
#mpirun ---mca btl openib,self -np 160 ./hw.exe



