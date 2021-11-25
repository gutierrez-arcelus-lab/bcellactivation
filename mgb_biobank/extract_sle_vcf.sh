#!/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

LABSHARE=/lab-share/IM-Gutierrez-e2/Public/ 
BED=${LABSHARE}/vitor/ase/mgb_biobank/sle_variants/sle_hg38.bed 

for CHR in {1..22}
    do
	for BATCH in {01..03}
	do
	    VCFIN=$LABSHARE/mgb_biobank/04${BATCH}/vcf/chr${CHR}.dose.vcf.gz 
	    VCFOUT=$TEMP_WORK/VCF/chr${CHR}.04${BATCH}.sle.vcf.gz

	    bcftools view -R $BED $VCFIN |\
		bcftools annotate -x INFO,FORMAT -O z -o $VCFOUT -

	    tabix -p vcf $VCFOUT
	done

	VCFCHR=$TEMP_WORK/VCF/chr${CHR}.04{01..03}.sle.vcf.gz 
	MERGEDCHR=$TEMP_WORK/VCF/chr${CHR}.sle.merged.vcf.gz
	bcftools merge -O z -o $MERGEDCHR $( eval echo $VCFCHR )
	tabix -p vcf $MERGEDCHR
done

ALLVCFS=${TEMP_WORK}/VCF/chr{1..22}.sle.merged.vcf.gz 
MERGED=${LABSHARE}/vitor/ase/mgb_biobank/sle_variants/sle.MGB.vcf
bcftools concat $( eval echo $ALLVCFS ) |\
    bcftools norm -m +both -O v -o $MERGED -

#rm ${TEMP_WORK}/VCF/chr*.sle*.vcf.gz*
