#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=24gb
#SBATCH --time=24:00:00 
#SBATCH -p bch-compute 
#SBATCH --array=1-216
#SBATCH --job-name=mbv
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/mbv-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/mbv-%A-%a.err

source /programs/biogrids.shrc
export SAMTOOLS_X=1.13
export QTLTOOLS_X=1.3.1

# INPUT
SPEC=${SLURM_SUBMIT_DIR}/array_spec_mbv.tsv
PREFIX=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SPEC )
BATCH=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $SPEC )

BAM=${TEMP_WORK}/bam/highinput_231227/${PREFIX}_Aligned.sortedByCoord.out.bam
VCF=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/bcell_rnaseq/old/3-mbv/data/allchr.${BATCH}.vcf.gz

samtools index -@ $SLURM_CPUS_PER_TASK $BAM 

OUTDIR=${SLURM_SUBMIT_DIR}/results
mkdir -p $OUTDIR

QTLtools mbv \
    --bam $BAM \
    --vcf $VCF \
    --out ${OUTDIR}/${PREFIX}_${BATCH}.txt \
    --filter-mapping-quality 255 \
    --filter-minimal-coverage 20 \
    --filter-imputation-qual 0 \
    --filter-imputation-prob 0 
