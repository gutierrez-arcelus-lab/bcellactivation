#!/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

VCFIN=( $(echo ./VCF/chr{1..22}.1000G.vcf.gz) )
VCFIN+=( $(echo ./VCF/chrX.1000G.vcf.gz ) )
ALLVCF=$(echo ${VCFIN[*]})
OUT=./VCF/allchr.1000G.vcf.gz

bcftools concat -o $OUT -O z $ALLVCF
tabix -p vcf $OUT
