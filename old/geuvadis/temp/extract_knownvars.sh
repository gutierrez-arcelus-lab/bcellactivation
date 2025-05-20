#! /usr/bin/env bash

source /programs/biogrids.shrc

SAMPLE=ERR188022
SAMPLE_1KG=$( awk -v var="$SAMPLE" '$0 ~ var {print $1}' ../E-GEUV-1.sdrf.txt | uniq )
DIR=../results/gatk
VCF=${DIR}/${SAMPLE}.filtered.vcf.gz

## 1000 Genomes remapped to GRCh38
#for i in {1..22}
#do
#    CHR=$i
#    VCF1KG=../../../vcf_1000G/ALL.chr${CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
#    VCF1KG_FILT=${DIR}/${SAMPLE}_chr${CHR}.1000KG.vcf.gz
#    POS=${DIR}/${SAMPLE}.chr${CHR}.pos
#
#    zcat $VCF |\
#	grep -v "^#" |\
#	awk '{ print $1"\t"$2 }' |\
#	awk 'gsub(/chr/, "", $0)' | \
#	awk -v chr="$CHR" '$1 == chr' > $POS
#
#    bcftools view -R $POS $VCF1KG |\
#	bcftools view -s $SAMPLE_1KG - |\
#	bcftools annotate -x INFO -O z -o $VCF1KG_FILT 
#done
#
#CHRVCFS=${DIR}/${SAMPLE}_chr{1..22}.1000KG.vcf.gz 
#
#bcftools concat -o ${DIR}/${SAMPLE}_1000KG.vcf.gz $(eval echo $CHRVCFS) 
#
#rm $(eval echo $CHRVCFS)
#rm ${DIR}/${SAMPLE}.chr*.pos

# NYGC 1000 Genomes
#for i in {1..22}
#do
#    CHR=chr$i
#    VCF1KG=/reference_databases/1000G_VCF/GRCh38/Phased_VCFs/CCDG_14151_B01_GRM_WGS_2020-08-05_${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz
#    VCF1KG_FILT=${TEMP_WORK}/${SAMPLE}.${CHR}.nygc1000g.vcf.gz
#    POS=${TEMP_WORK}/${SAMPLE}.${CHR}.pos
#
#    zcat $VCF |\
#	grep -v "^#" |\
#	awk '{ print $1"\t"$2 }' |\
#	awk -v chr="$CHR" '$1 == chr' > $POS
#
#    bcftools view -R $POS $VCF1KG |\
#	bcftools view -s $SAMPLE_1KG - |\
#	bcftools annotate -x INFO -O z -o $VCF1KG_FILT 
#done

CHRVCFS=${TEMP_WORK}/${SAMPLE}.chr{1..22}.nygc1000g.vcf.gz 

bcftools concat -o ../results/gatk/${SAMPLE}.nygc1000g.vcf.gz $(eval echo $CHRVCFS) 

