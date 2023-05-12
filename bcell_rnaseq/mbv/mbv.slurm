#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=16gb
#SBATCH --time=24:00:00 
#SBATCH -p bch-compute 
#SBATCH --array=102
#SBATCH --job-name=MBV
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/MBV-%A-%a

source /programs/biogrids.shrc
export SAMTOOLS_X=1.13
export QTLTOOLS_X=1.3.1

# INPUT
META=${SLURM_SUBMIT_DIR}/mbv.spec
PREFIX=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $META )
BATCH=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $META )

BAM=${TEMP_WORK}/${PREFIX}.uniq.mkdups.wasp.bam 
VCF=${SLURM_SUBMIT_DIR}/data/allchr.snps.maf.${BATCH}.vcf.gz

#samtools index -@ $SLURM_CPUS_PER_TASK $BAM 

QTLtools mbv \
    --bam $BAM \
    --vcf $VCF \
    --out ${SLURM_SUBMIT_DIR}/results/${PREFIX}_${BATCH}.txt \
    --filter-mapping-quality 255 \
    --filter-base-quality 0 \
    --filter-minimal-coverage 10 \
    --filter-imputation-qual 0 \
    --filter-imputation-prob 0 
