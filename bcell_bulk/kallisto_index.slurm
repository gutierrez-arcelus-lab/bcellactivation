#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=kallisto-index
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /home/ch229163/ase/log/kallisto-index.%j

source /programs/biogrids.shrc
export KALLISTO_X=0.46.2

GENCODE=$HOME/ase/data/gencode.transcripts.fa
INDEX=$HOME/ase/indices/kallisto

mkdir -p $INDEX

kallisto index -i $INDEX/transcripts.idx $GENCODE

