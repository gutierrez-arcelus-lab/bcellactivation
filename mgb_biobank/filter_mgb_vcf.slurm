#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-198
#SBATCH --job-name=vcf_process
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/vcf_process.%A.%a

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Define directory and chromosome given array design file
DESIGN=${SLURM_SUBMIT_DIR}/array_design.txt
BATCH=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $DESIGN )
CHR=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $DESIGN )

# Define IO
LABSHARE=/lab-share/IM-Gutierrez-e2/Public
VCFIN=${LABSHARE}/mgb_biobank/$BATCH/vcf/chr${CHR}.dose.vcf.gz 
IDS=${SLURM_SUBMIT_DIR}/results/females_${BATCH}.txt

OUTDIR=${SLURM_SUBMIT_DIR}/results/VCF
VCFOUT=${OUTDIR}/chr${CHR}.MGB.${BATCH}.vcf.gz

# Select biallelic SNPs only, with MAF >= 1%;
# remove variants with any missing genotypes;
# extract GT info and discard dosages and probabilities;
bcftools view --samples-file $IDS --force-samples $VCFIN |\
    bcftools annotate -x INFO,FORMAT - |\
    bcftools norm -m +both - |\
    bcftools view --genotype ^miss --types snps -m2 -M2 --min-af 0.01:minor - |\
    bgzip -c > $VCFOUT    

tabix -p vcf $VCFOUT
