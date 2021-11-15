#!/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Merge VCFs for each chromosome into a single VCF
OUTDIR=/lab-share/IM-Gutierrez-e2/Public/vitor/ase/mgb_biobank/results 
VCF=$OUTDIR/allchr.vcf.gz

bcftools concat -o $VCF -O z $(ls -v $TEMP_WORK/VCF/chr*.merged.pruned.vcf)
tabix -p vcf $VCF

# Make PLINK binary files
plink --vcf $VCF \
    --keep-allele-order \
    --make-bed \
    --out $OUTDIR/allchr

# Subset VCF for individuals in 1000G or MGB
VCF1000G=$OUTDIR/allchr.1000G.vcf.gz
VCFMGB=$OUTDIR/allchr.MGB.vcf.gz

grep "^NA[0-9]\|^HG[0-9]" $OUTDIR/allchr.fam |\
    awk '{ print $1"_"$2 }' - > $OUTDIR/samples.1000G.txt

grep "^[0-9]" $OUTDIR/allchr.fam |\
    awk '{ print $1"_"$2 }' - > $OUTDIR/samples.MGB.txt

bcftools view --samples-file $OUTDIR/samples.1000G.txt -O z -o $VCF1000G $VCF
bcftools view --samples-file $OUTDIR/samples.MGB.txt -O z -o $VCFMGB $VCF

plink --vcf $VCF1000G \
    --keep-allele-order \
    --make-bed \
    --out $OUTDIR/allchr.1000G

plink --vcf $VCFMGB \
    --keep-allele-order \
    --make-bed \
    --out $OUTDIR/allchr.MGB

