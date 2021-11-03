#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
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
OUTDIR=$TEMP_WORK/VCF

CHR=$SLURM_ARRAY_TASK_ID
VCF=$MGB/chr${CHR}.dose.vcf.gz
VCFOUT=$OUTDIR/chr${CHR}.filtered.vcf.gz
VARS=$OUTDIR/chr${CHR}.variants.txt.gz

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
DIRKG=$LABSHARE/vitor/vcf_1000G
VCFKG=$DIRKG/ALL.chr${CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
VCFKGOUT=$OUTDIR/chr${CHR}.1000G.filtered.vcf.gz
VARSKG=$OUTDIR/chr${CHR}.1000G.variants.txt.gz

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
POSFILE=$OUTDIR/chr${CHR}.filtered.pos
Rscript $HOME/ase/mgb_biobank/match_snps.R $VARS $VARSKG $POSFILE

# Filter MGB and 1000G VCFs for the same set of variants
# and merge files
VCFOUT2=$OUTDIR/chr${CHR}.filtered.v2.vcf.gz
VCFKGOUT2=$OUTDIR/chr${CHR}.1000G.filtered.v2.vcf.gz
VCFMERGED=$OUTDIR/chr${CHR}.merged.vcf.gz

tabix -p vcf $VCFOUT
tabix -p vcf $VCFKGOUT
bcftools view -R $POSFILE $VCFOUT -O z -o $VCFOUT2 
bcftools view -R $POSFILE $VCFKGOUT -O z -o $VCFKGOUT2

tabix -p vcf $VCFOUT2
tabix -p vcf $VCFKGOUT2
bcftools merge $VCFOUT2 $VCFKGOUT2 |\
    bcftools annotate -x 'INFO' - -O z -o $VCFMERGED

GDS=$OUTDIR/chr${CHR}.gds
Rscript $HOME/ase/mgb_biobank/prune_vcf.R $VCFMERGED $GDS
