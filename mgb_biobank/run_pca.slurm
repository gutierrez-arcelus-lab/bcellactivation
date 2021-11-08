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

# Merge VCFs for each chromosome into a single VCF
OUTDIR=/lab-share/IM-Gutierrez-e2/Public/vitor/ase/mgb_biobank/results 
VCF=$OUTDIR/allchr.vcf.gz

#bcftools concat -o $VCF -O z $(ls -v $TEMP_WORK/VCF/chr*.merged.pruned.vcf.gz)
#tabix -p vcf $VCF

# PCA
OUT=$OUTDIR/mgb

plink --vcf $VCF --pca --out $OUT
