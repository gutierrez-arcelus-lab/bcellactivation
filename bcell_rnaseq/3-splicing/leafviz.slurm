#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=4gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH --array=1-6
#SBATCH --job-name=leafviz
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/leafviz.%A.%a.out
#SBATCH -e /temp_work/ch229163/log/leafviz.%A.%a.err

source /programs/biogrids/biogrids.shrc
export PYTHON_X=2.7.2
export R_X=4.1

cd $SLURM_SUBMIT_DIR

CONTRAST=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' ./data/contrasts.txt )
COUNT=./data/clusters_perind_numers.counts.gz
GROUP=./data/groups_${CONTRAST}.tsv

/lab-share/IM-Gutierrez-e2/Public/tools/leafcutter/leafviz/prepare_results.R \
    --code ${CONTRAST} \
    --meta_data_file $GROUP \
    $COUNT \
    ./results/${CONTRAST}_cluster_significance.txt \
    ./results/${CONTRAST}_effect_sizes.txt \
    ./data/leafviz_annot \
    -o ./results/${CONTRAST}.Rdata
