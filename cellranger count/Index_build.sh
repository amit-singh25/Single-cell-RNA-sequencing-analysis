#!/bin/bash
#PBS -N index_build
#PBS -j oe 
##PBS -l file=32GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8	
#PBS -o ./ndex_build_$PBS_JOBID.out
#PBS -e ./ndex_build_$PBS_JOBID.err
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


#####downlaod 

wget ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz


wget ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz
gunzip Mus_musculus.GRCm38.93.gtf.gz

gn=~/genome/10x_genome

cellranger mkgtf ${gn}/Mus_musculus.GRCm38.93.gtf Mus_musculus.GRCm38.93.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense \
                 --attribute=gene_biotype:IG_LV_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:IG_V_pseudogene \
                 --attribute=gene_biotype:IG_D_gene \
                 --attribute=gene_biotype:IG_J_gene \
                 --attribute=gene_biotype:IG_J_pseudogene \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_C_pseudogene \
                 --attribute=gene_biotype:TR_V_gene \
                 --attribute=gene_biotype:TR_V_pseudogene \
                 --attribute=gene_biotype:TR_D_gene \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_J_pseudogene \
                 --attribute=gene_biotype:TR_C_gene


cellranger mkref --genome=mm10 \
                 --fasta=${gn}/Mus_musculus.GRCm38.dna.primary_assembly.fa \
                 --genes=${gn}/Mus_musculus.GRCm38.93.filtered.gtf \
                 --ref-version=3.0.0
