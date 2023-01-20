#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=PCA
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/PCA_expression.%j


# BED file with expression data
PHENO=$HOME/lupus/phenotypes.bed.gz

# Results
OUT=$HOME/lupus/data/pca
mkdir -p $OUT

# Run QTLtools pca on expression data
QTLtools pca --bed $PHENO --out $OUT/pheno_pcs --center --scale

