#!/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

LABSHARE=/lab-share/IM-Gutierrez-e2/Public/ 
BED=$LABSHARE/vitor/ase/mgb_biobank/sle_variants/sle_hg38.bed 
MERGED=$LABSHARE/vitor/ase/mgb_biobank/sle_variants/sle.MGB.vcf

for BATCH in {01..03}
    do
	for CHR in {1..22}
	do
	    VCFIN=$LABSHARE/mgb_biobank/04${BATCH}/vcf/chr${CHR}.dose.vcf.gz 
	    VCFOUT=$TEMP_WORK/VCF/chr${CHR}.04${BATCH}.sle.MGB.vcf.gz
	    VCFOUTLIST+=($VCFOUT)

	    bcftools view -R $BED $VCFIN |\
		bcftools annotate -x INFO,FORMAT -O z -o $VCFOUT -
	done
done

bcftools concat -O v -o $MERGED "${VCFOUTLIST[@]}"
