#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH --job-name=bcfconcat
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/bcftools_concat_%j.out
#SBATCH -e /temp_work/ch229163/log/bcftools_concat_%j.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

VCFS=( $(echo ${TEMP_WORK}/vcf/tmp/chr{{1..22},X}_merged.vcf.gz ) )
VCFOUT=${SLURM_SUBMIT_DIR}/data/allchroms_mgb.vcf.gz

bcftools concat $VCFS |\
    bcftools annotate -x QUAL,FILTER,INFO,FORMAT -O z -o $VCFOUT -

tabix -p vcf $VCFOUT

# Create VCF with variants only (no individual information)
bcftools view -G -O z -o ${SLURM_SUBMIT_DIR}/data/allchroms_mgb_variants.vcf.gz $VCFOUT
