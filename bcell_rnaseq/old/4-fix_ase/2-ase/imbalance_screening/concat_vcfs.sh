#!/usr/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

DIR=${TEMP_WORK}/vcf/langefeld
VCFS=( $(eval echo ${DIR}/chr{1..22}.MGB.merged.vcf.gz) )
VCFOUT=./data/allchr.mgb.vcf.gz 

bcftools concat "${VCFS[@]}" |\
    bcftools annotate -x FILTER,INFO,FORMAT -O z -o $VCFOUT -

tabix -f -p vcf $VCFOUT
