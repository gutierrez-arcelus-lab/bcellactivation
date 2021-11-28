#!/bin/bash

#SBATCH -N 1 
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=72:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-10
#SBATCH --job-name=ADMIXTURE-projection
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/ADMIXTURE-projection.%A.%a

K=$SLURM_ARRAY_TASK_ID

cd ${SLURM_SUBMIT_DIR}/results

cp allchr.merged.pruned.1000G.${K}.P allchr.merged.pruned.MGB.${K}.P.in 

admixture --cv -P allchr.merged.pruned.MGB.bed $K | tee admix.MGBprojection.${K}.log
