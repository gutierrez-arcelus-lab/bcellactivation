#!/usr/bin/bash

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Define IO
DB=/lab-share/IM-Gutierrez-e2/Public/References/dbSNP/GCF_000001405.39.gz 
OUT=./data/gwas_rsid_dbsnp.vcf

# CMD
bcftools view --threads 4 -r NC_000012.12 $DB |\
    bcftools annotate --threads 4 -x INFO - |\
    bcftools view --threads 4 --include 'ID=@data/rsids_bentham.txt' -O v -o $OUT -

