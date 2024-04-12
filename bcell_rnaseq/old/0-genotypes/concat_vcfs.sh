#!/usr/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

DIR=${TEMP_WORK}/vcf/ase 
VCFS=( $(eval echo ${DIR}/chr{{1..22},X}.MGB.merged.vcf.gz) )
VCFOUT=./data/allchr.mgb.vcf.gz 

bcftools concat -O z -o $VCFOUT "${VCFS[@]}"
tabix -f -p vcf $VCFOUT

META=./data/metadata.tsv

for i in {1..16}
do
    ID=$( awk 'FNR > 1' $META | awk -v ROW="$i" 'NR == ROW { print $2 }' )
    MGBID=$( awk 'FNR > 1' $META | awk -v ROW="$i" 'NR == ROW { print $3 }' )
    OUTI=./data/${ID}.vcf.gz 

    bcftools view --samples $MGBID --force-samples $VCFOUT |\
	bcftools view --genotype het -O z -o $OUTI -

    tabix -f -p vcf $OUTI 
done
