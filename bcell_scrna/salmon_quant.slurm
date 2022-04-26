#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=Salmon-quant
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/Salmon-quant.%j

source /programs/biogrids.shrc
export SALMON_X=1.5.1

# Alignment index
INDEX=$HOME/lupus/indices/SALMON_HUMAN_EBVDECOY

# Fastq file
FQ=${SLURM_SUBMIT_DIR}/mapping/unmapped_pilot_scrna_noTags.fastq.gz

# Results
OUT=${SLURM_SUBMIT_DIR}/mapping/salmon
mkdir -p $OUT

# Run Salmon quantification
salmon quant -i $INDEX -l A -r $FQ -p $SLURM_CPUS_PER_TASK -o $OUT \
    --gcBias --writeUnmappedNames
