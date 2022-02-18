#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-5
#SBATCH --job-name=Salmon-quant
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/Salmon-quant.%j

source /programs/biogrids.shrc
export SALMON_X=1.5.1

# Fastq files
LABSHARE=/lab-share/IM-Gutierrez-e2/Public 
FQDIR=$LABSHARE/B_cells/bulkTCpilot_1/34.198.31.178/210618_MG8989_fastq

# Alignment index
INDEX=$HOME/lupus/indices/SALMON

# Defining the samples
SAMPLELIST=/home/ch229163/lupus/bcell_bulk/bcell_samples.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )

# Defining the fastq files for each sample
FQ1=${FQDIR}/${SAMPLE}_R1_001.fastq.gz 
FQ2=${FQDIR}/${SAMPLE}_R2_001.fastq.gz 

# Results
OUT=$HOME/lupus/bcell_bulk/results/salmon/${SAMPLE}
mkdir -p $OUT

# Run Salmon quantification
salmon quant -i $INDEX -l A -1 $FQ1 -2 $FQ2 -p $SLURM_CPUS_PER_TASK -o $OUT \
    --seqBias --posBias --gcBias --validateMappings
