#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=SalmonIx
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/Salmon-index.%j

source /programs/biogrids.shrc
export SALMON_X=1.5.1

LABSHR=/lab-share/IM-Gutierrez-e2/Public
GENCODE=${LABSHR}/vitor/ase/data/gencode_plusEBV.transcripts.fa
INDEX=${LABSHR}/vitor/ase/indices

mkdir -p $INDEX

salmon index -t $GENCODE -i $INDEX/SALMON_HUMAN_EBV

