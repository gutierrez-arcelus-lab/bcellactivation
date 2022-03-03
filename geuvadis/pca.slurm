#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=PCA
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/PCA-geuvadis.%j

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12
export QTLTOOLS_X=1.3.1

PREFIX=${SLURM_SUBMIT_DIR}/VCF/allchr.1000G
VCF=${PREFIX}.vcf.gz

QTLtools pca \
    --vcf $VCF \
    --out $PREFIX \
    --center --scale --distance 5000
