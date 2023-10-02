#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=24gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -J Salmon-index
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/Salmon-index.%j.out
#SBATCH -e /temp_work/ch229163/log/Salmon-index.%j.err

source /programs/biogrids.shrc
export SALMON_X=1.10.1

GENCODE=/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.primary_assembly.transcripts.fa
INDEX=$HOME/lupus/indices

mkdir -p $INDEX

salmon index -t $GENCODE -i $INDEX/SALMON --tmpdir ${TEMP_WORK}/temp/salmon

