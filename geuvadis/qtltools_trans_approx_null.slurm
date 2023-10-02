#! /usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=QTLTRANS-permute
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/QTLTRANS-permute.%j

source /programs/biogrids.shrc
export QTLTOOLS_X=1.3.1
export R_X=4.1

GENOT=${SLURM_SUBMIT_DIR}/VCF/allchr.1000G.vcf.gz
#PHENO=${SLURM_SUBMIT_DIR}/ebv_phenotypes_corrected.bed.gz
PHENO=${SLURM_SUBMIT_DIR}/ebv_phenotypes_corrected_v2.bed.gz
OUTDIR=${SLURM_SUBMIT_DIR}/results/qtltools 
SEED=123
#OUTPERMUT=${OUTDIR}/trans.sample
OUTPERMUT=${OUTDIR}/trans.sample_v2

# Permutation pass
QTLtools trans \
    --vcf $GENOT \
    --bed $PHENO \
    --sample 86 \
    --out $OUTPERMUT

