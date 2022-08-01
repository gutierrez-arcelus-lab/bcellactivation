#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute-pe
#SBATCH --array=0-4
#SBATCH --job-name=leafcut
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/leafcutter.%A.%a

source /programs/biogrids/biogrids.shrc
export PYTHON_X=2.7.2
export R_X=4.1

cd $SLURM_SUBMIT_DIR

CELLTYPES=("aN" "DN" "rN" "SM" "T3")
TYPE="${CELLTYPES[$SLURM_ARRAY_TASK_ID]}"
OUT=./results/$TYPE
JUNC=${OUT}/juncfiles.txt

# Intron clustering
python $HOME/leafcutter/clustering/leafcutter_cluster_regtools.py \
    -j $JUNC \
    -m 30 \
    -l 100000 \
    -p 0.05 \
    -o ${OUT}/leafcutter

for f in $(cat juncfiles.txt) 
do 
    mv $(basename $f).leafcutter.sorted.gz ${OUT}/
done

# Differential intron excision analysis
COUNT=${OUT}/leafcutter_perind_numers.counts.gz
GROUP=${OUT}/groups_file.txt
EXON=./grch38_exon.txt.gz
OUT=${OUT}/leafcutter

$HOME/leafcutter/scripts/leafcutter_ds.R \
    -p $SLURM_CPUS_PER_TASK \
    -e $EXON \
    -o $OUT \
    $COUNT \
    $GROUP
