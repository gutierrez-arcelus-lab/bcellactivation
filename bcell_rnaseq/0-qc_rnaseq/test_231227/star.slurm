#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=48gb
#SBATCH --time=24:00:00 
#SBATCH -p bch-compute,bch-compute-pe 
#SBATCH --array=2-25
#SBATCH --job-name=STAR
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/STAR-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/STAR-%A-%a.err

source /programs/biogrids.shrc
export STAR_X=2.7.11a

# INPUT
META=${SLURM_SUBMIT_DIR}/metadata.tsv
DONORID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $META )
REPID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $META )
STIM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $3 }' $META )

PREFIX=${DONORID}_${REPID}_${STIM} 

FQ1=${TEMP_WORK}/fastq/highinput_231227/${PREFIX}_val_1.fq.gz
FQ2=${TEMP_WORK}/fastq/highinput_231227/${PREFIX}_val_2.fq.gz

# INDEX
INDEX=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/indices/STAR

# OUTPUT
OUT=${TEMP_WORK}/bam/highinput_231227/${PREFIX} 

# Run mapping
STAR --runMode alignReads \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $INDEX \
    --readFilesIn $FQ1 $FQ2 \
    --readFilesCommand zcat  \
    --outFilterMismatchNoverReadLmax 0.04 \
    --twopassMode Basic \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${OUT}_
