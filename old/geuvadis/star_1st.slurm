#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=48gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-2
#SBATCH --job-name=STAR-1st
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/STAR-1st-geuvadis.o%A_%a

source /programs/biogrids.shrc
export STAR_X=2.7.9a

SAMPLELIST=$HOME/ase/geuvadis/samples.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )
FQDIR=/reference_databases/1000G_Geuvadis_rnaseq_lymphoblastoid_cell_lines/fastq_files

FQ1=$FQDIR/${SAMPLE}_1.fastq.gz
FQ2=$FQDIR/${SAMPLE}_2.fastq.gz

LABSHARE=/lab-share/IM-Gutierrez-e2/Public/vitor 
INDEX=$LABSHARE/ase/indices/starindex
OUTDIR=$LABSHARE/ase/geuvadis/results/star/pass_1st
OUT=$OUTDIR/${SAMPLE}

mkdir -p $OUTDIR

STAR --runMode alignReads \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir $INDEX \
    --readFilesIn $FQ1 $FQ2 \
    --readFilesCommand zcat  \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix ${OUT}_

