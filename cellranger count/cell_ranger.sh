#!/bin/bash
#PBS -N cell_ranger
#PBS -j oe 
##PBS -l file=32GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8	
#PBS -o ./cell_ranger_$PBS_JOBID.out
#PBS -e ./cell_ranger_$PBS_JOBID.err
echo $PBS_JOBID
echo $PBS_JOBNAME
cd $PBS_O_WORKDIR

#echo "=========================================================="
#echo "Starting on : $(date)"
#echo "Running on node : $(hostname)"
#echo "Current directory : $(pwd)"
#echo "Current job ID : $JOB_ID"
#echo "Current job name : $JOB_NAME"
#echo "Task index number : $SGE_TASK_ID"
#echo "=========================================================="
