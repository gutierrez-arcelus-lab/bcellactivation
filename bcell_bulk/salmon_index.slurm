#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=Salmon-index
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /home/ch229163/ase/log/Salmon-index.%j

source /programs/biogrids.shrc
export SALMON_X=1.5.1

GENCODE=$HOME/ase/data/gencode.transcripts.fa
INDEX=$HOME/ase/indices

mkdir -p $INDEX

salmon index -t $GENCODE -i $INDEX/SALMON

