#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH -p bch-largemem
#SBATCH --array=34,38%1
#SBATCH --job-name=coloc
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/coloc-lupus-%a

source /programs/biogrids.shrc
export R_X=4.2.0

cd $SLURM_SUBMIT_DIR

DATASET=bentham
REGIONID=$SLURM_ARRAY_TASK_ID 

Rscript run_coloc_tier2.R $DATASET $REGIONID
  
