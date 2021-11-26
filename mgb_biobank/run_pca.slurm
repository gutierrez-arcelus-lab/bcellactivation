#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=PCA
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/pca.%j

source /programs/biogrids.shrc

DIR=/lab-share/IM-Gutierrez-e2/Public/vitor/ase/mgb_biobank/results 
PREFIX=allchr.merged.pruned
VCF=${DIR}/${PREFIX}.vcf.gz
OUT=${DIR}/${PREFIX}.pca

plink --vcf $VCF --pca --out $OUT
