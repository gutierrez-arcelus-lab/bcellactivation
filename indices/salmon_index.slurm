#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=Salmon-index
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/Salmon-index.%j

source /programs/biogrids.shrc
export SALMON_X=1.5.1

GENCODE=$HOME/lupus/data/gencode.transcripts.fa
INDEX=$HOME/lupus/indices

mkdir -p $INDEX

salmon index -t $GENCODE -i $INDEX/SALMON

