#!/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

MERGEDCHR=( $(echo ./results/VCF/chr{1..22}.sle.merged.vcf.gz) )
MERGEDCHR+=( $(echo ./results/VCF/chrX.sle.merged.vcf.gz) )
ALLVCF=$(echo "${MERGEDCHR[*]}")

MERGED=./sle_variants/sle.MGB.vcf

bcftools concat $ALLVCF |\
    bcftools norm -m +both |\
    bcftools view -v snps -O v -o $MERGED -
