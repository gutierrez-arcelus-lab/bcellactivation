#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=48gb
#SBATCH --time=48:00:00 
#SBATCH -p bch-compute 
#SBATCH --array=1-74
#SBATCH --job-name=STAR
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/STAR-%A-%a.out
#SBATCH -e logs/STAR-%A-%a.err

source /programs/biogrids.shrc
export STAR_X=2.7.11a

cd $SLURM_SUBMIT_DIR

# INPUT
TARGET_ROW=$((SLURM_ARRAY_TASK_ID + 1))
META="../0-qc/data/metadata.tsv"

DONORID=$( awk -v row="$TARGET_ROW" 'FNR==row { print $1 }' $META )
REPID=$( awk -v row="$TARGET_ROW" 'FNR==row { print $2 }' $META )
STIM=$( awk -v row="$TARGET_ROW" 'FNR==row { print $3 }' $META )

PREFIX="${DONORID}_${REPID}_${STIM}" 

FQ1="${TEMP_WORK}/fastq/highinput_merged/${PREFIX}_R1_val_1.fq.gz"
FQ2="${TEMP_WORK}/fastq/highinput_merged/${PREFIX}_R2_val_2.fq.gz"

# INDEX
INDEX="./index"

# OUTPUT
OUT="${TEMP_WORK}/bam/highinput_merged/${PREFIX}" 

# Run mapping
STAR --runMode alignReads \
    --runThreadN "$SLURM_CPUS_PER_TASK" \
    --genomeDir "$INDEX" \
    --readFilesIn "$FQ1" "$FQ2" \
    --readFilesCommand zcat  \
    --outFilterMismatchNoverReadLmax 0.04 \
    --twopassMode Basic \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${OUT}_"
