#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH -p bch-compute
#SBATCH --job-name=ADMIX-projection
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /home/ch229163/log/ADMIX-projection.%o

K=5

cd $SLURM_SUBMIT_DIR/results

cp allchr.refpanel.$K.P allchr.MGB.$K.P.in 

admixture -P -j4 allchr.MGB.bed $K | tee admix.MGBprojection.$K.log
