#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=12gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH --array=1-9 
#SBATCH --job-name=leafcut
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/leafcutter_%A_%a.out
#SBATCH -o /temp_work/ch229163/log/leafcutter_%A_%a.err

source /programs/biogrids/biogrids.shrc
export PYTHON_X=2.7.2
export R_X=4.1

TOOLS=/lab-share/IM-Gutierrez-e2/Public/tools

# Input/Output
INFO=${SLURM_SUBMIT_DIR}/slurm_info.txt
DATASET=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $INFO )
CELLTYPE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $INFO )

JUNC=${SLURM_SUBMIT_DIR}/results/${DATASET}/${CELLTYPE}/juncfiles.txt
OUT=${SLURM_SUBMIT_DIR}/results/${DATASET}/${CELLTYPE}

# Intron clustering
python ${TOOLS}/leafcutter/clustering/leafcutter_cluster_regtools.py \
    -j $JUNC \
    -m 30 \
    -l 500000 \
    -p 0.05 \
    -r ${OUT} \
    -o leafcutter

rm ${OUT}/*.leafcutter.sorted.gz

# Differential splicing analysis for each dataset and cell type
COUNT=${OUT}/leafcutter_perind_numers.counts.gz
GROUP=${OUT}/groups_file.txt
EXON=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/splicing/grch38_exon.txt.gz

${TOOLS}/leafcutter/scripts/leafcutter_ds.R \
    -i 3 -g 3 \
    -p $SLURM_CPUS_PER_TASK \
    -e $EXON \
    -o ${OUT}/leafcutter \
    $COUNT \
    $GROUP

# Prepare results for visualization with leafviz
${TOOLS}/leafcutter/leafviz/prepare_results.R \
    --code leafcutter-${CELLTYPE} \
    --meta_data_file ${OUT}/groups_file.txt \
    $COUNT \
    ${OUT}/leafcutter_cluster_significance.txt \
    ${OUT}/leafcutter_effect_sizes.txt \
    ../leafviz/gencode_hg38/gencode_hg38 \
    -o ${OUT}/${CELLTYPE}.Rdata
