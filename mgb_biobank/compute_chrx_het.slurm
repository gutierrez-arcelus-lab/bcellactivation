#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=0-8
#SBATCH --job-name=chrXhet
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/chrX-het-genotyped.%A.%a

# Load tools
source /programs/biogrids.shrc
export VCFTOOLS_X=0.1.12b
export BCFTOOLS_X=1.12

# Define IO
BATCHES=( $(eval echo 04{01..08} 0410) )
BATCH="${BATCHES[$SLURM_ARRAY_TASK_ID]}"

LABSHARE=/lab-share/IM-Gutierrez-e2/Public
VCFIN=${LABSHARE}/mgb_biobank/${BATCH}/vcf/chrX.dose.vcf.gz
VCFGENOT=${TEMP_WORK}/VCF/chrX.${BATCH}.genotyped.vcf.gz
OUT=${TEMP_WORK}/VCF/chrX.${BATCH}.genotyped

if [ "$SLURM_ARRAY_TASK_ID" -lt 8 ]; then

    INFO=${LABSHARE}/mgb_biobank/${BATCH}/vcf/chrX.info.gz
    GENOTX=${TEMP_WORK}/VCF/chrX.${BATCH}.genotyped.txt
    zcat $INFO | awk 'FNR > 1 && $8 == "Genotyped" { print $1 }' > $GENOTX

elif [ "$SLURM_ARRAY_TASK_ID" -eq 8 ]; then

    GENOTX=${LABSHARE}/mgb_biobank/${BATCH}/vcf/chrX.genotyped.txt
fi

bcftools view --include ID==@$GENOTX -O z -o $VCFGENOT $VCFIN

vcftools --gzvcf $VCFGENOT --het --out $OUT
