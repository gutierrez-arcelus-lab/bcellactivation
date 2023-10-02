#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-74
#SBATCH --job-name=FASTQC
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/FASTQC.%A-%a.log
#SBATCH -e /temp_work/ch229163/log/FASTQC.%A-%a.err

source /programs/biogrids.shrc
export CUTADAPT_X=4.1
export FASTQC_X=0.11.9
export PYTHON_X=3.8.8

# INPUT
DIR=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_rnaseq 
META=${DIR}/metadata.tsv
FQDIR=${DIR}/34.198.31.178/221025_MG10430_fastq
FQ1=${FQDIR}/$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $4 }' )
FQ2=${FQDIR}/$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $5 }' )

# OUTPUT
OUT=${TEMP_WORK}/fastq/highinput
mkdir -p $OUT

trim_galore --paired \
    --length 50 \
    --fastqc \
    --cores $SLURM_CPUS_PER_TASK \
    -o $OUT \
    $FQ1 $FQ2
