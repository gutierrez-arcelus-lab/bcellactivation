#! /usr/bin/env bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# 1000 Genomes SNPs
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx

# dbSNP from GATK
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# Newer dbSNP from NCBI
wget https://ftp.ncbi.nlm.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz
wget https://ftp.ncbi.nlm.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz.tbi
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt 

Rscript prepare_gatk_files.R

bcftools annotate --rename-chrs ./dbsnp_names.txt -o dbsnp_155.hg38.vcf.gz -O z GCF_000001405.39.gz

# Indels from GATK
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

