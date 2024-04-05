#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-110
#SBATCH --job-name=bcftools
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/subset-mgb-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/subset-mgb-%A-%a.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

SPEC=${SLURM_SUBMIT_DIR}/array_spec.txt
BATCH=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $1 }' $SPEC )
CHR=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' $SPEC )

VCFIN=/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank/${BATCH}/vcf/${CHR}.dose.vcf.gz 
IDS=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/bcell_rnaseq/4-fix_ase/0-genotypes/data/${BATCH}.txt
VARS=${SLURM_SUBMIT_DIR}/data/gwas_vars.tsv
VCFOUT=${TEMP_WORK}/vcf/langefeld/${CHR}.MGB.${BATCH}.vcf.gz

bcftools view --threads $SLURM_CPUS_PER_TASK -R $VARS $VCFIN |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --samples-file $IDS --force-samples - |\
    bcftools norm --threads $SLURM_CPUS_PER_TASK -m +both -O z -o $VCFOUT -

tabix -f -p vcf $VCFOUT
