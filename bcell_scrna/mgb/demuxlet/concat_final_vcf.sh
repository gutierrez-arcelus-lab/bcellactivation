#!/usr/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

VCFS=${TEMP_WORK}/chr{1..22}.MGB.merged.vcf.gz
VCFOUT=./data/allchr.mgb.vcf.gz 

bcftools concat -o $VCFOUT -O z $( eval echo $VCFS )
tabix -f -p vcf $VCFOUT
