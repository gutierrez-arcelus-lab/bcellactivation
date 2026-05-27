#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH --array=1-92
#SBATCH --job-name=bcftools
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/bcftools_mgb_%A_%a.out
#SBATCH -e /temp_work/ch229163/log/bcftools_mgb_%A_%a.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

SPEC=${SLURM_SUBMIT_DIR}/array_spec.txt
CHROM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SPEC )
BATCH=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $SPEC )
IDS=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $3 }' $SPEC )

VCFIN=/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/${BATCH}/vcf/chr${CHROM}.dose.vcf.gz 
VCFOUT=${TEMP_WORK}/vcf/tmp/chr${CHROM}_${BATCH}.vcf.gz

bcftools view --threads $SLURM_CPUS_PER_TASK --samples $IDS $VCFIN |\
    bcftools view --threads $SLURM_CPUS_PER_TASK -i 'INFO/R2>=0.9' -O z -o $VCFOUT -

tabix -p vcf $VCFOUT
