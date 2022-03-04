#! /usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=QTLTRANS
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/QTLTRANS-adjust.%j

source /programs/biogrids.shrc
export QTLTOOLS_X=1.3.1

GENOT=${SLURM_SUBMIT_DIR}/VCF/allchr.1000G.vcf.gz
PHENO=${SLURM_SUBMIT_DIR}/ebv_phenotypes_corrected.bed.gz
OUTDIR=${SLURM_SUBMIT_DIR}/results/qtltools 
NOMOUT=${OUTDIR}/trans.adjust

QTLtools trans \
    --vcf $GENOT \
    --bed $PHENO \
    --adjust ${OUTDIR}/trans.sample.best.txt.gz \
    --threshold 0.99 \
    --out $NOMOUT

