#! /usr/bin/env bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12
export R_X=4.1

DBSNP_IN=$TEMP_WORK/GCF_000001405.39.gz
DBSNP_OUT=$TEMP_WORK/dbsnp_155.hg38.vcf.gz 

# Newer dbSNP from NCBI
wget -P $TEMP_WORK https://ftp.ncbi.nlm.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz
wget -P $TEMP_WORK https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt 

# convert chr names to Gencode standard
Rscript prepare_gatk_files.R $TEMP_WORK/GCF_000001405.39_GRCh38.p13_assembly_report.txt  

bcftools annotate \
    --rename-chrs $TEMP_WORK/chr_dbsnpToGencode_names.txt \
    -O z -o $DBSNP_OUT \
    $DBSNP_IN

tabix -p vcf $DBSNP_OUT

# select chrs in GRCh38 Primary Assembly
DBSNP_PRI=./dbsnp_155.hg38pri.vcf.gz 
bcftools view -R ref_pri_chr.bed -O z -o $DBSNP_PRI $DBSNP_OUT 
tabix -p vcf $DBSNP_PRI

#KGP_IN=$TEMP_WORK/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
#KGP_OUT=./1000g.gatkbundle.hg38.vcf.gz
#
#INDEL_IN=$TEMP_WORK/Homo_sapiens_assembly38.known_indels.vcf.gz
#INDEL_OUT=$TEMP_WORK/knownindels_gatkbundle.vcf.gz
#
#MILLS_IN=$TEMP_WORK/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
#MILLS_OUT=$TEMP_WORK/mills_1000g_indels_gatkbundle.vcf.gz
#
## 1000 Genomes SNPs
#wget -P $TEMP_WORK https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
#
#bgzip -c $KGP_IN > $KGP_OUT
#tabix -p vcf $KGP_OUT
#
## 1000 Genomes file from GATK is corrupted;
## File interrupted on chr15, wrong chr lengths, etc
## There are forum posts discussing those issues
## Using original from 1000G Project instead:
#
#wget -P $TEMP_WORK http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz
#
## Indels from GATK
#wget -P $TEMP_WORK https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
#wget -P $TEMP_WORK https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
#
## dbSNP from GATK
#wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
#wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
#

#
#bcftools annotate \
#    --rename-chrs $TEMP_WORK/chr_UcscToGencode_names.txt \
#    -O z -o $INDEL_OUT \
#    $INDEL_IN
#
#bcftools annotate \
#    --rename-chrs $TEMP_WORK/chr_UcscToGencode_names.txt \
#    -O z -o $MILLS_OUT \
#    $MILLS_IN
#
#tabix -p vcf $INDEL_OUT
#tabix -p vcf $MILLS_OUT
#
#INDEL_PRI=./knownindels_gatkbundle.hg38pri.vcf.gz 
#MILLS_PRI=./mills_1000g_indels_gatkbundle.hg38pri.vcf.gz 
#
#
#bcftools view -R ref_pri_chr.bed -O z -o $INDEL_PRI $INDEL_OUT 
#tabix -p vcf $INDEL_PRI
#
#bcftools view -R ref_pri_chr.bed -O z -o $MILLS_PRI $MILLS_OUT 
#tabix -p vcf $MILLS_PRI
#
