#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=32gb
#SBATCH --time=24:00:00 
#SBATCH -p bch-compute 
#SBATCH --array=1-666
#SBATCH --job-name=mbv
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/mbv-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/mbv-%A-%a.err

source /programs/biogrids.shrc
export SAMTOOLS_X=1.13
export QTLTOOLS_X=1.3.1

# INPUT
SPEC=${SLURM_SUBMIT_DIR}/array_spec_mbv.txt
PREFIX=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SPEC )
BATCH=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $SPEC )

BAM=${TEMP_WORK}/bam/highinput/${PREFIX}.uniq.mkdups.wasp.bam
VCF=${SLURM_SUBMIT_DIR}/data/allchr.${BATCH}.vcf.gz

samtools index -@ $SLURM_CPUS_PER_TASK $BAM 

OUTDIR=${SLURM_SUBMIT_DIR}/results
mkdir -p $OUTDIR

QTLtools mbv \
    --bam $BAM \
    --vcf $VCF \
    --out ${OUTDIR}/${PREFIX}_${BATCH}.txt \
    --filter-mapping-quality 255 \
    --filter-base-quality 0 \
    --filter-minimal-coverage 20 \
    --filter-imputation-qual 0 \
    --filter-imputation-prob 0 
