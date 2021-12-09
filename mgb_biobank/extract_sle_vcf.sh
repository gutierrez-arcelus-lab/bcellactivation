#!/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# convert Hg19 coordinations to Hg38
BED=./sle_variants/sle.bed 
CHAIN=/reference_databases/ReferenceGenome/liftover_chain/hg19/hg19ToHg38.over.chain.gz
OUT=./sle_variants/sle_hg38.bed
FAIL=./sle_variants/sle_failtolift.bed

liftOver $BED $CHAIN $OUT $FAIL 

# Extract SLE variants from MGB VCFs
LABSHARE=/lab-share/IM-Gutierrez-e2/Public/ 
BED=${LABSHARE}/vitor/ase/mgb_biobank/sle_variants/sle_hg38.bed 

for CHR in $( seq 1 22; echo X )
    do
	for BATCH in {01..07}
	do
	    VCFIN=${LABSHARE}/mgb_biobank/04${BATCH}/vcf/chr${CHR}.dose.vcf.gz 
	    VCFOUT=${TEMP_WORK}/VCF/chr${CHR}.04${BATCH}.sle.vcf.gz
	    SAMPLES=./results/females_04${BATCH}.txt

	    bcftools view -R $BED --samples-file $SAMPLES --force-samples $VCFIN |\
		bcftools annotate -x INFO,FORMAT -O z -o $VCFOUT -

	    tabix -p vcf $VCFOUT
	done

	VCFCHR=${TEMP_WORK}/VCF/chr${CHR}.04{01..07}.sle.vcf.gz 
	MERGEDCHR=${TEMP_WORK}/VCF/chr${CHR}.sle.merged.vcf.gz
	MERGEDLIST+=($MERGEDCHR)

	bcftools merge -O z -o $MERGEDCHR $( eval echo $VCFCHR )
	tabix -p vcf $MERGEDCHR
done

MERGED=${LABSHARE}/vitor/ase/mgb_biobank/sle_variants/sle.MGB.vcf

bcftools concat "${MERGEDLIST[@]}" |\
    bcftools norm -m +both -O v -o $MERGED -

#rm ${TEMP_WORK}/VCF/chr*.sle*.vcf.gz*

