#!/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

LABSHARE=/lab-share/IM-Gutierrez-e2/Public/ 
DIR=$LABSHARE/mgb_biobank/0401/vcf
BED=$LABSHARE/vitor/ase/mgb_biobank/sle_variants/sle_hg38.bed 
OUT=$LABSHARE/vitor/ase/mgb_biobank/sle_variants/sle.MGB.vcf

for i in {1..22}
do
    bcftools view -R $BED $DIR/chr$i.dose.vcf.gz |\
	bcftools annotate -x INFO,FORMAT -O z -o $TEMP_WORK/VCF/chr$i.sle.MGB.vcf.gz - 
done

bcftools concat -O v -o $OUT $(ls -v $TEMP_WORK/VCF/chr{1..22}.sle.MGB.vcf.gz)
