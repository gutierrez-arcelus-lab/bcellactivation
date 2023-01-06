#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --ntasks=4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1
#SBATCH --job-name=coloc
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/coloc-lupus-%a

source /programs/biogrids.shrc
export R_X=4.2.0

REGION=$SLURM_ARRAY_TASK_ID
CPUS=$SLURM_CPUS_PER_TASK

cd $SLURM_SUBMIT_DIR

Rscript run_coloc.R $REGION $CPUS
