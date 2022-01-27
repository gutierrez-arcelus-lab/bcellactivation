#!/usr/bin/env bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

VCF=../results/gatk/ERR188022.filtered.vcf.gz
OUT=./ERR188022.norm.vcf.gz

bcftools norm -m -both $VCF |\
    bcftools norm -m +both -O z -o $OUT -
