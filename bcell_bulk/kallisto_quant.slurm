#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-5
#SBATCH --job-name=kallisto-quant
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /home/ch229163/ase/log/kallisto-quant.%j

source /programs/biogrids.shrc
export KALLISTO_X=0.46.2

LABSHARE=/lab-share/IM-Gutierrez-e2/Public 
FQDIR=$LABSHARE/B_cells/bulkTCpilot_1/34.198.31.178/210618_MG8989_fastq

INDEX=$HOME/ase/indices/kallisto/transcripts.idx
SAMPLELIST=/home/ch229163/ase/bcell_bulk/bcell_samples.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )

FQ1=${FQDIR}/${SAMPLE}_R1_001.fastq.gz 
FQ2=${FQDIR}/${SAMPLE}_R2_001.fastq.gz 

OUT=$HOME/ase/bcell_bullk/results/kallisto/${SAMPLE}

mkdir -p $OUT

kallisto quant -i $INDEX -o $OUT -b 100 -t $SLURM_CPUS_PER_TASK --bias $FQ1 $FQ2
