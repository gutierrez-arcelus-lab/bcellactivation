#!/usr/bin/bash

source /programs/biogrids/biogrids.shrc
export BCFTOOLS_X=1.12

VCF19=/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.25.gz
REG19=./pos_23andme.txt
OUT19=./dbsnp_23andme.vcf

VCF38=/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.39.gz
REG38=./pos_23andme_hg38.txt
OUT38=./dbsnp_23andme_hg38.vcf

#bcftools view --threads 4 \
#    -R $REG19 \
#    --types snps \
#    -O v -o $OUT19 $VCF19

bcftools view --threads 4 \
    -R $REG38 \
    --types snps \
    -O v -o $OUT38 $VCF38
