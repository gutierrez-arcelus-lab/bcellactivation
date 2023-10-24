#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=0-22
#SBATCH --job-name=vcf-merge
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/vcf-merge-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/vcf-merge-%A-%a.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Params
CHRS=( $(seq 1 22; echo X) )
CHR=chr${CHRS[$SLURM_ARRAY_TASK_ID]}
CPUS=$SLURM_CPUS_PER_TASK

DIR=${TEMP_WORK}/vcf/quasar

VCFS=( $(eval echo ${DIR}/${CHR}_04{03..07}.vcf.gz) )
VCFMERGED=${DIR}/${CHR}.merged.vcf.gz

#bcftools merge --threads $CPUS "${VCFS[@]}" |\
#    bcftools +fill-tags -O z -o $VCFMERGED - 

bcftools view -G -O z -o ${DIR}/${CHR}.merged.info.vcf.gz $VCFMERGED
