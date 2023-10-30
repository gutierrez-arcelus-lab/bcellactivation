#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=24gb
#SBATCH --time=24:00:00 
#SBATCH -p bch-compute 
#SBATCH --array=1-74
#SBATCH --job-name=QuASAR
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/QuASAR-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/QuASAR-%A-%a.err

source /programs/biogrids/biogrids.shrc
export SAMTOOLS_X=1.13
export BEDTOOLS_X=2.31.0
export R_X=4.1

LABSHR=/lab-share/IM-Gutierrez-e2/Public
PILEUPtoQUASAR=${LABSHR}/vitor/software/QuASAR/scripts/convertPileupToQuasar.R

META=${LABSHR}/vitor/lupus/bcell_rnaseq/4-fix_ase/1-mapping/metadata.tsv
SAMPLEID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' $META )
STIM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $3 }' $META )
PREFIX=${SAMPLEID}_${STIM}
BAM=${TEMP_WORK}/bam/fixase/${PREFIX}.uniq.mkdups.wasp.bam

GENOME=${LABSHR}/References/Genomes/hsapiens/GRCh38/GRCh38.primary_assembly.genome.fa
SNPS=${SLURM_SUBMIT_DIR}/data/snps_af.bed

samtools mpileup -f $GENOME -l $SNPS $BAM |\
    awk -v OFS='\t' '{ if ($4>0 && $5 !~ /[^\^][<>]/ && $5 !~ /\+[0-9]+[ACGTNacgtn]+/ && $5 !~ /-[0-9]+[ACGTNacgtn]+/ && $5 !~ /[^\^]\*/) print $1,$2-1,$2,$3,$4,$5,$6}' |\
    sortBed -i stdin |\
    intersectBed -a stdin -b $SNPS -wo |\
    cut -f 1-7,11-14 |\
    gzip > ${TEMP_WORK}/quasar/${PREFIX}.pileup.bed.gz 

MINCOV=10

Rscript $PILEUPtoQUASAR ${TEMP_WORK}/quasar/${PREFIX}.pileup.bed.gz $MINCOV

# QuASAR script writes output to SUBMIT_DIR
# Move files to TEMP dir instead
mv ${SLURM_SUBMIT_DIR}/${PREFIX}.quasar.in.gz ${TEMP_WORK}/quasar/ 
rm ${TEMP_WORK}/quasar/${PREFIX}.pileup.bed.gz 
