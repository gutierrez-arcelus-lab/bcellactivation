#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=0-5
#SBATCH --job-name=vcftools-het
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/vcftools-het.%A.%a

# Load tools
source /programs/biogrids.shrc
export VCFTOOLS_X=0.1.12b

# Define IO
BATCHES=( $(eval echo 04{01..06}) )
BATCH="${BATCHES[$SLURM_ARRAY_TASK_ID]}"

LABSHARE=/lab-share/IM-Gutierrez-e2/Public
VCFIN=${LABSHARE}/mgb_biobank/${BATCH}/vcf/chrX.dose.vcf.gz

OUT=${TEMP_WORK}/VCF/chrX.MGB.${BATCH}

vcftools --gzvcf $VCFIN --het --out $OUT
