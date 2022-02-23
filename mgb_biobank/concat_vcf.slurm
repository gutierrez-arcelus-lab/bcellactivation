#!/bin/bash

#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=concat_vcf
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/concat_vcf.%j

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Merge VCFs for each chromosome into a single VCF
OUTDIR=${SLURM_SUBMIT_DIR}/results 
VCFIN=${TEMP_WORK}/VCF/chr{1..22}.merged.pruned.vcf
PREFIX=${OUTDIR}/allchr.merged.pruned 
VCFOUT=${PREFIX}.vcf.gz

bcftools concat -o $VCFOUT -O z $( eval echo $VCFIN )
tabix -p vcf $VCFOUT

# Make PLINK binary files
plink --vcf $VCFOUT \
    --keep-allele-order \
    --make-bed \
    --out $PREFIX

## Subset VCF for individuals in 1000G or MGB
## I used to do this in order to run ADMIXTURE
#VCF1000G=${PREFIX}.1000G.vcf.gz
#VCFMGB=${PREFIX}.MGB.vcf.gz
#
#grep "^NA[0-9]\|^HG[0-9]" ${PREFIX}.fam |\
#    awk '{ print $1"_"$2 }' - > ${OUTDIR}/samples.1000G.txt
#
#grep "^[0-9]" ${PREFIX}.fam |\
#    awk '{ print $1"_"$2 }' - > ${OUTDIR}/samples.MGB.txt
#
#bcftools view --samples-file ${OUTDIR}/samples.1000G.txt -O z -o $VCF1000G $VCFOUT
#bcftools view --samples-file ${OUTDIR}/samples.MGB.txt -O z -o $VCFMGB $VCFOUT
#
#plink --vcf $VCF1000G \
#    --keep-allele-order \
#    --make-bed \
#    --out ${PREFIX}.1000G
#
#plink --vcf $VCFMGB \
#    --keep-allele-order \
#    --make-bed \
#    --out ${PREFIX}.MGB
#
