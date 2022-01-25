#! /usr/bin/env bash

SAMPLE=ERR188022
SAMPLE_1KG=$( awk -v var="$SAMPLE" '$0 ~ var {print $1}' E-GEUV-1.sdrf.txt | uniq )
DIR=./results/gatk
VCF=${DIR}/${SAMPLE}.filtered.denorm.vcf.gz

# 1000 Genomes remapped to GRCh38
for i in {1..22}
do
    CHR=$i
    VCF1KG=../../../vcf_1000G/ALL.chr${CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
    VCF1KG_FILT=${DIR}/${SAMPLE}_chr${CHR}.1000KG.vcf.gz
    POS=${DIR}/${SAMPLE}.chr${CHR}.pos

    zcat $VCF |\
	grep -v "^#" |\
	awk '{ print $1"\t"$2 }' |\
	awk 'gsub(/chr/, "", $0)' | \
	awk -v chr="$CHR" '$1 == chr' > $POS

    bcftools view -R $POS $VCF1KG |\
	bcftools view -s $SAMPLE_1KG - |\
	bcftools annotate -x INFO -O z -o $VCF1KG_FILT 
done

CHRVCFS=${DIR}/${SAMPLE}_chr{1..22}.1000KG.vcf.gz 

bcftools concat -o ${DIR}/${SAMPLE}_1000KG.vcf.gz $(eval echo $CHRVCFS) 

rm $(eval echo $CHRVCFS)
rm ${DIR}/${SAMPLE}.chr*.pos

# DBSNP
DBSNP=../data/dbsnp_155.hg38pri.vcf.gz
ALLPOS=${DIR}/${SAMPLE}.pos
DBSNP_OUT=${DIR}/${SAMPLE}_dbsnp.vcf.gz

zcat $VCF |\
    grep -v "^#" |\
    awk '{ print $1"\t"$2 }' |\
    uniq > $ALLPOS

bcftools view -R $ALLPOS -O z -o $DBSNP_OUT $DBSNP

rm $ALLPOS
