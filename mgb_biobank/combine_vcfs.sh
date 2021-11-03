#!/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

VCFOUT=/lab-share/IM-Gutierrez-e2/Public/vitor/mgb_biobank/allchr.vcf.gz

bcftools concat -o $VCFOUT -O z $(ls -v $TEMP_WORK/VCF/chr*.pruned.vcf)
tabix -p vcf $VCFOUT
