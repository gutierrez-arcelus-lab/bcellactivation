#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-92
#SBATCH --job-name=SubsetVCF
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/subsetVCF-%A-%a

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

DESIGN=${SLURM_SUBMIT_DIR}/array_subset.tsv
BATCH=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $1 }' $DESIGN )
CHR=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' $DESIGN )

VCFIN=/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/${BATCH}/vcf/${CHR}.dose.vcf.gz 
IDS=${SLURM_SUBMIT_DIR}/data/${BATCH}.txt
VCFOUT=${TEMP_WORK}/vcf/mbv/${CHR}.MGB.${BATCH}.vcf.gz

bcftools view --threads $SLURM_CPUS_PER_TASK --samples-file $IDS --force-samples $VCFIN |\
    bcftools norm --threads $SLURM_CPUS_PER_TASK -m +both - |\
    bcftools view --threads $SLURM_CPUS_PER_TASK -i 'INFO/R2>=0.8' --genotype ^miss --types snps -m2 -M2 -O z -o $VCFOUT -

tabix -f -p vcf $VCFOUT
