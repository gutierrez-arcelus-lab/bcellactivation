#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=concat_vcf
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/concat_vcf.%j

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Merge VCFs for each chromosome into a single VCF
OUTDIR=${SLURM_SUBMIT_DIR}/results/VCF
VCFIN=${OUTDIR}/chr{1..22}.isec.merged.pruned.vcf
PREFIX=${OUTDIR}/allchr.isec.merged.pruned 
VCFOUT=${PREFIX}.vcf.gz

bcftools concat --threads $SLURM_CPUS_PER_TASK -o $VCFOUT -O z $( eval echo $VCFIN )
tabix -p vcf $VCFOUT

# Make PLINK binary files
plink --vcf $VCFOUT \
    --keep-allele-order \
    --make-bed \
    --out $PREFIX

