#! /usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=32gb
#SBATCH --time=96:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=QTLTRANS
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/QTLTRANS.%j

source /programs/biogrids.shrc
export QTLTOOLS_X=1.3.1

GENOT=${SLURM_SUBMIT_DIR}/VCF/allchr.1000G.vcf.gz
PHENO=${SLURM_SUBMIT_DIR}/peer_residuals.bed
OUTDIR=${SLURM_SUBMIT_DIR}/results/qtltools 
NOMOUT=${OUTDIR}/trans.nominal

QTLtools trans \
    --vcf $GENOT \
    --bed $PHENO \
    --nominal --threshold 1e-5 \
    --out $NOMOUT

