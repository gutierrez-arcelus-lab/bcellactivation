#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=2gb
#SBATCH --time=6:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-115
#SBATCH --job-name=FASTQC
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/FASTQC.%A-%a.out
#SBATCH -e logs/FASTQC.%A-%a.err

source /programs/biogrids.shrc
export FASTQC_X=0.11.9

# INPUT
META=${SLURM_SUBMIT_DIR}/data/metadata.tsv
FQ1=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $4 }')
FQ2=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $5 }')

#OUTPUT
OUT=${TEMP_WORK}/fastq

trim_galore --paired \
    --length 75 \
    --fastqc \
    --cores $SLURM_CPUS_PER_TASK \
    -o $OUT \
    $FQ1 $FQ2
