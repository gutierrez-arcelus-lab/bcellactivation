#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=24gb
#SBATCH --time=24:00:00 
#SBATCH -p bch-compute 
#SBATCH --array=2-75
#SBATCH --job-name=Salmon
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/Salmon-%A-%a.out
#SBATCH -o /temp_work/ch229163/log/Salmon-%A-%a.err

source /programs/biogrids.shrc
export SALMON_X=1.5.1

# INPUT
META=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/bcell_rnaseq/0-qc_rnaseq/metadata.tsv
DONORID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $META )
REPID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $META )
STIM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $3 }' $META )
FQ1=${TEMP_WORK}/fastq/highinput_merged/${DONORID}_${REPID}_${STIM}_R1_val_1.fq.gz
FQ2=${TEMP_WORK}/fastq/highinput_merged/${DONORID}_${REPID}_${STIM}_R2_val_2.fq.gz

# INDEX
INDEX=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/indices/SALMON

# OUTPUT
OUT=${SLURM_SUBMIT_DIR}/results/${DONORID}_${REPID}_${STIM}
mkdir -p $OUT

# EXECUTE
salmon quant \
    -i $INDEX -l A \
    -1 $FQ1 -2 $FQ2 \
    -p $SLURM_CPUS_PER_TASK \
    -o $OUT \
    --seqBias --posBias --gcBias

