#! /usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=32gb
#SBATCH --time=96:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=PEER
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/PEER.%j

source /programs/biogrids.shrc
export R_X=4.1

cd $SLURM_SUBMIT_DIR

Rscript peer.R
