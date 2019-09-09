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


genome=~/genome/10x_genome/mouse_genome
data=~/single_cell/fastqs
name=$1

cellranger count --chemistry=SC3Pv3 \
--id=${name}-REX --project=P180721 \
--transcriptome=${genome} \
--fastqs=${data} \
--sample=${name} \
--localcores=7 \
--localmem=128



#cellranger count --chemistry=SC3Pv3 \
#--id=633818_28-REX --project=P180721 \
#--transcriptome=~/genome/10x_genome/mouse_genome \
#--fastqs=~/single_cell/fastqs \
#--sample=633818_28 \
#--localcores=7 \
#--localmem=128

cellranger count --id=control_2 --sample=63381822 --fastqs=~/single_cell/fastqs/ --transcriptome=~/genome/10x_genome/mouse_genome --expect-cells=500

