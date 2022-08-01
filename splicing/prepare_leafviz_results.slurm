#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute-pe
#SBATCH --array=0-4
#SBATCH --job-name=prepleaf
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/prep-leafviz.%A.%a

source /programs/biogrids.shrc
export R_X=4.1

cd $SLURM_SUBMIT_DIR

CELLTYPES=("aN" "DN" "rN" "SM" "T3")
TYPE="${CELLTYPES[$SLURM_ARRAY_TASK_ID]}"

DIR=./results/${TYPE}

$HOME/leafcutter/leafviz/prepare_results.R \
    --code leafcutter-${TYPE} \
    --meta_data_file ${DIR}/groups_file.txt \
    ${DIR}/leafcutter_perind_numers.counts.gz \
    ${DIR}/leafcutter_cluster_significance.txt \
    ${DIR}/leafcutter_effect_sizes.txt \
    ./leafviz/gencode_hg38/gencode_hg38 \
    -o ${DIR}/${TYPE}.Rdata
