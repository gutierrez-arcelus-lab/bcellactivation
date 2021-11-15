#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=ADMIXTURE-refpanel
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /home/ch229163/log/ADMIXTURE-refpanel.%o

DIR=/lab-share/IM-Gutierrez-e2/Public/vitor/ase/mgb_biobank/results
PREFIX=allchr.refpanel
K=5

admixture $DIR/$PREFIX.bed $K -j4 | tee $DIR/admix.refpanel.log

mv $SLURM_SUBMIT_DIR/$PREFIX.$K.P $SLURM_SUBMIT_DIR/$PREFIX.$K.Q $DIR/
