#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=48gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -J PLINKPCA
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/pca-%j.out
#SBATCH -e /temp_work/ch229163/log/pca-%j.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12
export PLINK2_X=2.00a3

# Merge VCFs for each chromosome into a single VCF
OUTDIR=${TEMP_WORK}/vcf/prs
VCFIN=${OUTDIR}/chr{1..22}.merged.pruned.vcf
PREFIX=${OUTDIR}/allchr.merged.pruned 
VCFOUT=${PREFIX}.vcf.gz
PCAOUT=${PREFIX}.plink.pca

bcftools concat --threads $SLURM_CPUS_PER_TASK -O z -o $VCFOUT $( eval echo $VCFIN )
tabix -p vcf $VCFOUT

# Run PCA
plink2 --vcf $VCFOUT --pca approx --threads $SLURM_CPUS_PER_TASK --out $PCAOUT
