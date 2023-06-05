#!/usr/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

VCFS=( $(eval echo ${TEMP_WORK}/vcf/mbv/chr{{1..22},X}.MGB.merged.vcf.gz) )
VCFOUT=./data/allchr.r2filtered.mgb.vcf.gz 

bcftools concat -O z -o $VCFOUT "${VCFS[@]}"
tabix -f -p vcf $VCFOUT

META=./data/metadata.tsv

for i in {1..5};
do
    ID=$( awk -v ROW="$i" 'NR == ROW { print $2 }' $META )
    MGBID=$( awk -v ROW="$i" 'NR == ROW { print $3 }' $META )
    OUTI=./data/${ID}.vcf.gz 
    
    bcftools view --samples $MGBID --force-samples $VCFOUT |
	bcftools view --genotype het -O z -o $OUTI -

    tabix -f -p vcf $OUTI 
done
