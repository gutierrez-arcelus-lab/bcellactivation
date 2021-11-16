#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --time=120:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-10
#SBATCH --job-name=ADMIXTURE-cv
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /home/ch229163/log/ADMIXTURE-cv.%A.%a

DIR=/lab-share/IM-Gutierrez-e2/Public/vitor/ase/mgb_biobank/results
K=$SLURM_ARRAY_TASK_ID

admixture --cv $DIR/allchr.1000G.bed $K -j4 | tee $DIR/admix.1000G.cv$K.log

mv $SLURM_SUBMIT_DIR/allchr.1000G.$K.P $SLURM_SUBMIT_DIR/allchr.1000G.$K.Q $DIR/
