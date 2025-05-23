#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=36gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-5
#SBATCH --job-name=STAR-1st
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/STAR-1st.o%A_%a

source /programs/biogrids.shrc
export STAR_X=2.7.9a

# Input
LABSHR=/lab-share/IM-Gutierrez-e2/Public
INDEX=${LABSHR}/vitor/lupus/indices/starindex
FQDIR=${LABSHR}/B_cells/bulkTCpilot_1/34.198.31.178/210618_MG8989_fastq
SAMPLELIST=${LABSHR}/vitor/lupus/bcell_bulk/bcell_samples.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )

FQ1=${FQDIR}/${SAMPLE}_R1_001.fastq.gz
FQ2=${FQDIR}/${SAMPLE}_R2_001.fastq.gz

# Output
OUTDIR=${TEMP_WORK}/star/pass_1st
OUT=${OUTDIR}/${SAMPLE}

mkdir -p $OUTDIR

# Command
STAR --runMode alignReads \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $INDEX \
    --readFilesIn $FQ1 $FQ2 \
    --readFilesCommand zcat  \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outFilterType BySJout \
    --outSJfilterReads Unique \
    --outSJfilterCountUniqueMin 6 2 2 2 \
    --outSJfilterCountTotalMin 6 2 2 2 \
    --outSAMtype None \
    --outFileNamePrefix ${OUT}_

