#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=4gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=0-22
#SBATCH --job-name=VCFinfo 
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/VCFinfo.%A.%a 

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Directories
LABSHR=/lab-share/IM-Gutierrez-e2/Public 
OUTDIR=${TEMP_WORK}

# Chromosome
ALLCHRS=( $(seq 1 22; echo X) )
CHR=chr${ALLCHRS[$SLURM_ARRAY_TASK_ID]}

#for BATCH in {0401..0408}; 
#do
#    INFO=${LABSHR}/mgb_biobank/${BATCH}/vcf/${CHR}.info.gz
#    OUT=${OUTDIR}/${CHR}_${BATCH}_typed.txt
#
#    zcat $INFO | awk 'FNR > 1 && $8 ~ /Genotyped/ {print $1}' > $OUT
#done

VCF10=${LABSHR}/mgb_biobank/0410/vcf/${CHR}.dose.vcf.gz
OUT10=${OUTDIR}/${CHR}_0410_typed.txt

zcat $VCF10 | awk '$0 !~ /^#/ && $8 ~ /TYPED/ {print $3}' > $OUT10
