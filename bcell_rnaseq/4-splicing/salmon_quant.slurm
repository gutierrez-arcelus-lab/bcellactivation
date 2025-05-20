#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH --array=1-74
#SBATCH --job-name=Salmon
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/Salmon-low.%A-%a.out
#SBATCH -e /temp_work/ch229163/log/Salmon-low.%A-%a.err

source /programs/biogrids.shrc
export SALMON_X=1.10.1

# INPUT
META=${SLURM_SUBMIT_DIR}/metadata_qced.tsv
SAMPLEID=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $1 }' )
FQ1=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' ) 
FQ2=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $3 }' ) 

FQLIST1=( $(echo $FQ1 | tr "," " ") )
FQLIST2=( $(echo $FQ2 | tr "," " ") )

# Index
INDEX=${SLURM_SUBMIT_DIR}/index/salmon

# Output
OUT=${SLURM_SUBMIT_DIR}/salmon_quant/${SAMPLEID}
mkdir -p $OUT

# Run Salmon quantification
salmon quant -i $INDEX -l A \
    -1 ${FQLIST1[@]} -2 ${FQLIST2[@]} \
    -p $SLURM_CPUS_PER_TASK \
    -o $OUT \
    --seqBias --posBias --gcBias
