#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --time=48:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-10
#SBATCH --job-name=ADMIXTURE-cv
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/ADMIXTURE-cv.%A.%a

PREFIX=allchr.merged.pruned.1000G 
K=$SLURM_ARRAY_TASK_ID

cd $SLURM_SUBMIT_DIR/results

admixture --cv ${PREFIX}.bed $K -j4 | tee ${PREFIX}.cv${K}.log

