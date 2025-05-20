#! /usr/bin/env bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=32gb
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
PHENO=${SLURM_SUBMIT_DIR}/ebv_phenotypes_corrected.bed.gz
OUTDIR=${SLURM_SUBMIT_DIR}/results/qtltools 
SEED=123
OUTPERMUT=${OUTDIR}/trans.perm.${SEED}

# Permutation pass
QTLtools trans \
    --vcf $GENOT \
    --bed $PHENO \
    --permute --threshold 1e-5 \
    --seed $SEED \
    --out $OUTPERMUT

# FDR
#NOMOUT=${OUTDIR}/trans.nominal
#OUT=${OUTDIR}/trans.fdr.txt
#
#Rscript ${SLURM_SUBMIT_DIR}/qtltools_runFDR_ftrans.R \
#    ${NOMOUT}.hits.txt.gz \
#    ${OUTPERMUT}.hits.txt.gz \
#    $OUT
