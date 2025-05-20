#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-8,10
#SBATCH -J chrx-hom
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/chrx-hom-%A-%a.out
#SBATCH -e /temp_work/ch229163/log/chrx-hom-%A-%a.err

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12
export VCFTOOLS_X=0.1.12b

# Define IO
LABSHARE=/lab-share/IM-Gutierrez-e2/Public
CHR=chrX
BATCH=04$(printf "%02d" $SLURM_ARRAY_TASK_ID)
VCFIN=${LABSHARE}/External_datasets/mgb_biobank/$BATCH/vcf/${CHR}.dose.vcf.gz 
IDS=${SLURM_SUBMIT_DIR}/mgb_data/ids_per_batch/${BATCH}.txt
OUTPREFIX=${TEMP_WORK}/vcf/prs/chrx/chrX.MGB.${BATCH}
VCFOUT=${OUTPREFIX}.vcf.gz

# CMD
bcftools view --threads $SLURM_CPUS_PER_TASK -r chrX:3000000-155000000 --samples-file $IDS --force-samples $VCFIN |\
    bcftools view --threads $SLURM_CPUS_PER_TASK -i 'INFO/R2>=0.9 && INFO/IMPUTED==0' - |\
    bcftools annotate --threads $SLURM_CPUS_PER_TASK -x FILTER,INFO,FORMAT -O z -o $VCFOUT -

tabix -p vcf $VCFOUT

# Compute chrX heterozygosity
vcftools --gzvcf $VCFOUT --het --out $OUTPREFIX
