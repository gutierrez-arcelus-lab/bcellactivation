#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute-pe
#SBATCH --array=0-5
#SBATCH --job-name=leafcutter
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/leafcutter.%A.%a

cd $SLURM_SUBMIT_DIR

CELLTYPES=( $(ls ./results) )
TYPE="${CELLTYPES[$SLURM_ARRAY_TASK_ID]}"

COUNT=./results/${TYPE}/${TYPE}_perind_numers.counts.gz
GROUP=./results/${TYPE}/groups_file.txt
EXON=./grch38_exon.txt.gz
OUT=./results/${TYPE}/leafcutter

$HOME/leafcutter/scripts/leafcutter_ds.R \
    -p $SLURM_CPUS_PER_TASK \
    -e $EXON \
    -o $OUT \
    $COUNT \
    $GROUP
