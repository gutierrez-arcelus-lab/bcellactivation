#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=01:00:00
#SBATCH -p bch-compute
#SBATCH --array=0-22
#SBATCH --job-name=VCFinfo 
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/VCFinfo.%A.%a 

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

LABSHARE=/lab-share/IM-Gutierrez-e2/Public 

# Chromosome
ALLCHRS=( $(seq 1 22; echo X) )
CHR=chr${ALLCHRS[$SLURM_ARRAY_TASK_ID]}
SLE=${SLURM_SUBMIT_DIR}/sle_variants/sle_variants_hg38.txt

VCFIN=$(echo ${LABSHARE}/mgb_biobank/0410/vcf/${CHR}.dose.vcf.gz)
OUTDIR=${SLURM_SUBMIT_DIR}/results/VCF
OUT=${OUTDIR}/${CHR}_0410_info.txt

bcftools view -H -R $SLE $VCFIN | awk '{print $3,$8}' OFS='\t' > $OUT
