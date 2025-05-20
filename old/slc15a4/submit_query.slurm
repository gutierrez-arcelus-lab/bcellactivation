#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-564
#SBATCH -J eqtl_cat
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/eqtl_cat-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/eqtl_cat-%A-%a.err

# Load tools
source /programs/biogrids.shrc
export R_X=4.1

cd ${SLURM_SUBMIT_DIR} 

META=./data/eqtl_datasets.tsv
DATASET=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' )

Rscript ./query_catalogue.R $DATASET
