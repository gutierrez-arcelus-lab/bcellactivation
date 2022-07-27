#!/usr/bin/bash

source /programs/biogrids/biogrids.shrc
export BCFTOOLS_X=1.12

VCF=/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.25.gz
OUT=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/CITEseq_pilot_2/dbsnp.vcf

bcftools view \
    --include ID==@/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/CITEseq_pilot_2/snps_23andMe.txt \
    -O v -o $OUT $VCF
