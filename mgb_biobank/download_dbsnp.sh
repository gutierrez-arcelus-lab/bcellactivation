#!/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

wget https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.39.gz
tabix -p vcf GCF_000001405.39.gz 

Rscript convert_dbsnp_names.R

POS=./sle_variants/sle_variants_hg38_dbSNP.txt
OUT=./sle_variants/sle_variants_hg38_dbSNP.vcf

bcftools view -R $POS GCF_000001405.39.gz -O z -o $OUT 
