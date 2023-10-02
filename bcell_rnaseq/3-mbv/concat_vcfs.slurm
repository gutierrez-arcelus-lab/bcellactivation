#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-8,10
#SBATCH --job-name=vcf-cat
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/vcf-concat-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/vcf-concat-%A-%a.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

BATCH=04$(printf "%02d" $SLURM_ARRAY_TASK_ID)
VCFS=( $(eval echo ${TEMP_WORK}/vcf/mbv/chr{{1..22},X}.MGB.${BATCH}.vcf.gz) )
VCFOUT=${SLURM_SUBMIT_DIR}/data/allchr.${BATCH}.vcf.gz

bcftools concat --threads $SLURM_CPUS_PER_TASK -O z -o $VCFOUT "${VCFS[@]}"
tabix -f -p vcf $VCFOUT
