#!/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# convert Hg19 coordinations to Hg38
LABSHARE=/lab-share/IM-Gutierrez-e2/Public 
BED=${LABSHARE}/vitor/ase/mgb_biobank/sle_variants/sle_langefeld_bentham.bed 
CHAIN=/reference_databases/ReferenceGenome/liftover_chain/hg19/hg19ToHg38.over.chain.gz

PREFIX="${BED%.*}"
BED38=${PREFIX}_hg38.bed
FAIL=${PREFIX}_failtolift.bed

liftOver $BED $CHAIN $BED38 $FAIL 

POS=${PREFIX}_hg38.txt
awk -F'\t' '{ printf ("%s\t%s\n", $1, $2) }' $BED38 > $POS

# Extract SLE variants from MGB VCFs
for CHR in $( seq 1 22; echo X )
    do
	for BATCH in {01..08}
	do
	    VCFIN=${LABSHARE}/mgb_biobank/04${BATCH}/vcf/chr${CHR}.dose.vcf.gz 
	    VCFOUT=${TEMP_WORK}/VCF/chr${CHR}.04${BATCH}.sle.vcf.gz
	    SAMPLES=./results/females_04${BATCH}.txt

	    bcftools view -R $POS --samples-file $SAMPLES --force-samples $VCFIN |\
		bcftools annotate -x INFO,FORMAT -O z -o $VCFOUT -

	    tabix -p vcf $VCFOUT

	    echo "chr ${CHR} batch ${BATCH} done!"
	done

	VCFCHR=${TEMP_WORK}/VCF/chr${CHR}.04{01..08}.sle.vcf.gz 
	MERGEDCHR=${TEMP_WORK}/VCF/chr${CHR}.sle.merged.vcf.gz
	MERGEDLIST+=($MERGEDCHR)

	echo "merging all batches for chr ${CHR}..."
	bcftools merge -O z -o $MERGEDCHR $( eval echo $VCFCHR )
	tabix -p vcf $MERGEDCHR
done

MERGED=${LABSHARE}/vitor/ase/mgb_biobank/sle_variants/sle.MGB.vcf

echo "Concatenating all VCFs..."
bcftools concat "${MERGEDLIST[@]}" |\
    bcftools norm -m +both -O v -o $MERGED -

#rm ${TEMP_WORK}/VCF/chr*.sle*.vcf.gz*

echo "Done!"
