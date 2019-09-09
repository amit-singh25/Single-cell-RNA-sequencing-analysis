#!/bin/bash
#PBS -N cellranger_aggr
##PBS -j oe
#PBS -m e
##PBS -l file=100GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   
#PBS -e ./cellranger_aggr_$PBS_JOBID.err           # stderr file
#PBS -o ./cellranger_aggr_$PBS_JOBID.out
#PBS -V
#echo $PBS_JOBNAME
#echo $PBS_JOBID

data=~/single_cell
name=$1

cellranger aggr --id=all_in_one --csv=${data}/agg_samples.csv --normalize=mapped
