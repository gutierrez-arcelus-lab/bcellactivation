#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=12gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-22
#SBATCH --job-name=vcf_process
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /home/ch229163/ase/log/vcf_process.%A.%a

# Load tools
source /programs/biogrids.shrc

export VCFTOOLS_X=0.1.12b
export BCFTOOLS_X=1.12
export R_X=4.1

# Specify input data
LABSHARE=/lab-share/IM-Gutierrez-e2/Public
MGB=$LABSHARE/mgb_biobank/0401/vcf

CHR=$SLURM_ARRAY_TASK_ID
VCF=$MGB/chr${CHR}.dose.vcf.gz
VCFOUT=$TEMP_WORK/chr${CHR}.filtered.vcf.gz
VARS=$TEMP_WORK/chr${CHR}.variants.txt.gz

# Select biallelic SNPs only, with MAF >= 1%;
# remove variants with any missing genotypes;
# extract GT info and discard dosages and probabilities;
# remove variants exceding p=0.05/79M for a HWE test.
bcftools view --genotype ^miss --types snps -m2 -M2 $VCF |\
    bcftools view --min-af 0.01:minor - |\
    bcftools annotate -x 'FORMAT' - |\
    vcftools --vcf - --hwe 0.00000000063 --recode --stdout |\
    bgzip -c > $VCFOUT

# Save Chromosome, Position, ID, Ref and Alt alleles to a file
zcat $VCFOUT |\
    grep -v '^##' |\
    awk '{print $1, $2, $3, $4, $5}' |\
    bgzip -c > $VARS

# Filter 1000G
DIRKG=/reference_databases/1000G_VCF/GRCh38/Genotype_VCFs
VCFKG=$DIRKG/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.recalibrated_variants.vcf.gz
VCFKGOUT=$TEMP_WORK/chr${CHR}.1000G.filtered.vcf.gz
VARSKG=$TEMP_WORK/chr${CHR}.1000G.variants.txt.gz

bcftools view --genotype ^miss --types snps -m2 -M2 $VCFKG |\
    bcftools view --min-af 0.01:minor - |\
    bcftools annotate -x 'FORMAT' - |\
    bgzip -c > $VCFKGOUT

# Save Chromosome, Position, ID, Ref and Alt alleles to a file
zcat $VCFKGOUT |\
    grep -v '^##' |\
    awk '{print $1, $2, $3, $4, $5}' |\
    bgzip -c > $VARSKG

# Match variants with those in 1000 Genomes data
Rscript $HOME/ase/mgb_biobank/match_snps.R $CHR

# Filter MGB and 1000G VCFs for the same set of variants
VCFOUT2=$TEMP_WORK/chr${CHR}.filtered.v2.vcf.gz
VCFKGOUT2=$TEMP_WORK/chr${CHR}.1000G.filtered.v2.vcf.gz

bcftools view -R $TEMP_WORK/chr${CHR}.filtered.pos $VCFOUT -O z -o $VCFOUT2 
bcftools view -R $TEMP_WORK/chr${CHR}.filtered.pos $VCFKGOUT -O z -o $VCFKGOUT2

#GDS=$TEMP_WORK/chr${CHR}.gds
#Rscript $HOME/ase/mgb_biobank/prune_vcf.R $VCF $GDS
