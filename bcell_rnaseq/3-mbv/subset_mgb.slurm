#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-207
#SBATCH --job-name=vcf-sub
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/vcf-subset-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/vcf-subset-%A-%a.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

SPEC=${SLURM_SUBMIT_DIR}/array_spec_vcf.txt
BATCH=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $1 }' $SPEC )
CHR=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' $SPEC )

VCFIN=/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/${BATCH}/vcf/${CHR}.dose.vcf.gz 
VCFOUT=${TEMP_WORK}/vcf/mbv/${CHR}.MGB.${BATCH}.vcf.gz

bcftools view --threads $SLURM_CPUS_PER_TASK -i "INFO/TYPED==1" --min-ac 1:minor $VCFIN |\
    bcftools norm --threads $SLURM_CPUS_PER_TASK -m +both - |\
    bcftools annotate --threads $SLURM_CPUS_PER_TASK -x FILTER,INFO,FORMAT -O z -o $VCFOUT -

tabix -f -p vcf $VCFOUT
