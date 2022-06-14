#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute-pe
#SBATCH --array=0-5
#SBATCH --job-name=DTU
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/DTU.%A.%a

export R_X=4.1

cd $SLURM_SUBMIT_DIR

CELLTYPES=( $(ls ./results) )
TYPE="${CELLTYPES[$SLURM_ARRAY_TASK_ID]}"

Rscript dtu.R $TYPE
