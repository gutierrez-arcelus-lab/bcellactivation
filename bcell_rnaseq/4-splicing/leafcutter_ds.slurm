#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH --array=0-2
#SBATCH --job-name=leafcut
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/leafcutter_ds.%A.%a.out
#SBATCH -e /temp_work/ch229163/log/leafcutter_ds.%A.%a.err

source /programs/biogrids/biogrids.shrc
export PYTHON_X=2.7.2
export R_X=4.1

cd $SLURM_SUBMIT_DIR

# IO
STIMS=( TLR7 BCR DN2 )
STIM=${STIMS[$SLURM_ARRAY_TASK_ID]}
COUNT=./data/clusters_perind_numers.counts.gz
GROUP=./data/groups_${STIM}.txt
EXON=./data/exon.txt.gz

/lab-share/IM-Gutierrez-e2/Public/tools/leafcutter/scripts/leafcutter_ds.R \
    -e $EXON \
    -p $SLURM_CPUS_PER_TASK \
    -o ./results/${STIM} \
    $COUNT \
    $GROUP
