#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=24gb
#SBATCH -t 08:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-4
#SBATCH -J homer
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/homer.%A-%a.out
#SBATCH -e /temp_work/ch229163/log/homer.%A-%a.err

source /programs/biogrids.shrc
export HOMER_X=4.11

cd ${SLURM_SUBMIT_DIR}

STIM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' ./stims.txt )

GENOME=/lab-share/IM-Gutierrez-e2/Public/References/Genomes/hsapiens/GRCh38/GRCh38.primary_assembly.genome.fa
BED=./data/${STIM}.bed
BG=./data/consensus.bed
OUT=./results/${STIM}

findMotifsGenome.pl $BED $GENOME $OUT \
    -bg $BG \
    -size 200 \
    -p $SLURM_CPUS_PER_TASK \
    -mask
