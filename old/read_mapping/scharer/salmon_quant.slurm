#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute-pe
#SBATCH --array=1-37
#SBATCH --job-name=SalmonQ
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/Salmon-quant.%j

source /programs/biogrids.shrc
export SALMON_X=1.5.1

# Fastq files
FQDIR=/lab-share/IM-Gutierrez-e2/Public/External_datasets/Scharer/RNAseq/fastq

# Alignment index
INDEX=$HOME/lupus/indices/SALMON_HUMAN_EBV

# Defining the samples
SAMPLELIST=${SLURM_SUBMIT_DIR}/samples_hc.txt
SAMPLE=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SAMPLELIST )

# Defining the fastq files for each sample
FQ1=${FQDIR}/${SAMPLE}_1.fastq.gz
FQ2=${FQDIR}/${SAMPLE}_2.fastq.gz 

# Results
OUT=${SLURM_SUBMIT_DIR}/salmon/${SAMPLE}
mkdir -p $OUT

# Run Salmon quantification
salmon quant -i $INDEX -l A -1 $FQ1 -2 $FQ2 -p $SLURM_CPUS_PER_TASK -o $OUT \
    --seqBias --posBias --gcBias
