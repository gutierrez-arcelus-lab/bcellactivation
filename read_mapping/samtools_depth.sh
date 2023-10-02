#!/usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH --array=1-227
#SBATCH --job-name=DEPTH
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/samtools-depth-%A-%a


source /programs/biogrids/biogrids.shrc
export SAMTOOLS_X=1.13

cd $SLURM_SUBMIT_DIR

LIST=./bam_list.txt

BAM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $LIST )
OUT=$(echo $BAM | sed 's|Aligned.sortedByCoord.out.bam|depth.txt|' )
BED=./eif4a2.bed

samtools depth -Q 255 -a -m 1000000 -b $BED $BAM > $OUT
