#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=24gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-174
#SBATCH -J Salmon
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/Salmon-quant.%A-%a.out
#SBATCH -e /temp_work/ch229163/log/Salmon-quant.%A-%a.err

source /programs/biogrids.shrc
export SALMON_X=1.10.1

# INPUT
META=${SLURM_SUBMIT_DIR}/data/metadata_pooledreps.tsv
SAMPLEID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $META)
FQ1=$(awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $META)
FQ2=$(awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $3 }' $META)
FQLIST1=(${FQ1//,/ })
FQLIST2=(${FQ2//,/ })
FQLIST1=( "${FQLIST1[@]/#/${TEMP_WORK}/fastq/lowinput/}" )
FQLIST2=( "${FQLIST2[@]/#/${TEMP_WORK}/fastq/lowinput/}" )

# Index
INDEX=${HOME}/lupus/indices/SALMON

# Output
OUT=${SLURM_SUBMIT_DIR}/results/salmon_pooledreps/${SAMPLEID}
mkdir -p $OUT

# Run Salmon quantification
salmon quant -i $INDEX -l A \
    -1 ${FQLIST1[@]} -2 ${FQLIST2[@]} \
    -p $SLURM_CPUS_PER_TASK \
    -o $OUT \
    --seqBias --posBias --gcBias \
    --numBootstraps 20
