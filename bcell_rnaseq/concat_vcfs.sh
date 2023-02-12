#!/usr/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

VCFS=( $(eval echo ${TEMP_WORK}/chr{1..22}.MGB.merged.vcf.gz) )
VCFS+=( ${TEMP_WORK}/chrX.MGB.merged.vcf.gz )
VCFOUT=./data/allchr.mgb.vcf.gz 

bcftools concat -O z -o $VCFOUT "${VCFS[@]}"
tabix -f -p vcf $VCFOUT

META=./data/metadata.tsv

for i in {1..16};
do
    ID=$( awk 'FNR > 1' $META | awk -v ROW="$i" 'NR == ROW { print $2 }' )
    MGBID=$( awk 'FNR > 1' $META | awk -v ROW="$i" 'NR == ROW { print $3 }' )
    
    bcftools view --samples $MGBID --force-samples $VCFOUT |
	bcftools view --genotype het -O z -o ./data/${ID}.vcf.gz -

    tabix -f -p vcf ./data/${ID}.vcf.gz 
done
