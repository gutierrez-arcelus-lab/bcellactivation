#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH --array=2-25
#SBATCH --job-name=fastqc
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/fastqc.%A-%a.log
#SBATCH -e /temp_work/ch229163/log/fastqc.%A-%a.err

source /programs/biogrids.shrc
export CUTADAPT_X=4.1
export FASTQC_X=0.11.9
export PYTHON_X=3.8.8

# INPUT
META=${SLURM_SUBMIT_DIR}/metadata.tsv
DONORID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $META )
REPID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $META )
STIM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $3 }' $META )
FQ1=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $4 }' $META )
FQ2=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $5 }' $META )

# OUTPUT
OUTDIR=${TEMP_WORK}/fastq/highinput_231227

# Run adapter trimming followed by FASTQC 
trim_galore --paired \
    --length 50 \
    --fastqc --fastqc_args "-t 4" \
    --cores $SLURM_CPUS_PER_TASK \
    -o $OUTDIR --basename ${DONORID}_${REPID}_${STIM} \
    $FQ1 $FQ2
