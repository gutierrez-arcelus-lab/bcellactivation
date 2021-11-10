#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=PCA
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /home/ch229163/ase/log/pca.%o

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

DIR=/lab-share/IM-Gutierrez-e2/Public/vitor/ase/mgb_biobank/results 
VCF=$DIR/allchr.vcf.gz
OUT=$DIR/mgb

plink --vcf $VCF --pca --out $OUT
