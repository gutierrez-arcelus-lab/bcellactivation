#!/usr/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.16

VCFS=${TEMP_WORK}/chr{1..22}.MGB.merged.vcf.gz
VCFOUT=../data/allchr.mgb.vcf.gz 

bcftools concat -O z -o $VCFOUT $( eval echo $VCFS )
tabix -f -p vcf $VCFOUT

# Sort chromosomes with the same order as BAM file
zcat $VCFOUT |\
    grep -v "^#" |\
    awk '{ print $1 }' |\
    sort |\
    uniq > ../data/chr_list.txt

VCFSORT=../data/allchr.mgb.sorted.vcf.gz 

cat ../data/chr_list.txt |\
    xargs tabix -h $VCFOUT |\
    bgzip > $VCFSORT

tabix -f -p vcf $VCFSORT

# Fix header, which does not reflect the sorting
Rscript fix_vcf_header.R

VCFSORTFIX=../data/allchr.mgb.sorted.reheader.vcf.gz 
bcftools reheader -h ../data/header.txt -o $VCFSORTFIX $VCFSORT
tabix -f -p vcf $VCFSORTFIX

rm ${VCFOUT}* ${VCFSORT}* 
rm ../data/chr_list.txt ../data/header.txt

# Extract gene regions
BED=../data/gene_regions.bed
OUT=../data/allchr.mgb.generegions.vcf.gz 

bcftools view -R $BED -O z -o $OUT $VCFSORTFIX 
tabix -f -p vcf $OUT
