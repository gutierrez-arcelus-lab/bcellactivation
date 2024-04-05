#!/usr/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

DIR=${TEMP_WORK}/vcf/susievars 
VCFS=( $(eval echo ${DIR}/chr{1,2,5,6,7,8,10,11,12,16,19,22}.MGB.merged.vcf.gz) )
VCFOUT=./allchr.mgb.vcf.gz 

bcftools concat -O z -o $VCFOUT "${VCFS[@]}"
tabix -f -p vcf $VCFOUT
