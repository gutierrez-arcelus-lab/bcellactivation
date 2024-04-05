#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=2-75
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
FQS1=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $4 }' $META )
FQS2=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $5 }' $META )

# OUTPUT
OUTDIR=${TEMP_WORK}/fastq/highinput_merged

# Merge fastqs
FQLIST1=(${FQS1//,/ })
FQLIST2=(${FQS2//,/ })

FQMERGE1=${OUTDIR}/${DONORID}_${REPID}_${STIM}_R1.fastq.gz
FQMERGE2=${OUTDIR}/${DONORID}_${REPID}_${STIM}_R2.fastq.gz

cat ${FQLIST1[@]} > $FQMERGE1
cat ${FQLIST2[@]} > $FQMERGE2

# Run adapter trimming followed by FASTQC 
trim_galore --paired \
    --length 50 \
    --fastqc \
    --cores $SLURM_CPUS_PER_TASK \
    -o $OUTDIR \
    $FQMERGE1 $FQMERGE2
