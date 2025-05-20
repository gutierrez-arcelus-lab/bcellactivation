#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -J bcf-concat
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/bcftools-concat-%j.out
#SBATCH -e /temp_work/ch229163/log/bcftools-concat-%j.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Merge VCFs for each chromosome into a single VCF
OUTDIR=${TEMP_WORK}/vcf/prs
VCFIN=${OUTDIR}/chr{{1..22},X}.merged.sentinels.vcf.gz
PREFIX=${OUTDIR}/allchr.mgb.sentinels 
VCFOUT=${PREFIX}.vcf.gz

bcftools concat --threads $SLURM_CPUS_PER_TASK -O z -o $VCFOUT $( eval echo $VCFIN )
tabix -p vcf $VCFOUT

