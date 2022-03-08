#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-22
#SBATCH --job-name=filter_1000G
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/filter_1000G.%A.%a

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Define IO
CHR=$SLURM_ARRAY_TASK_ID
LABSHARE=/lab-share/IM-Gutierrez-e2/Public
VCFIN=${LABSHARE}/vcf_1000G/ALL.chr${CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
OUTDIR=${SLURM_SUBMIT_DIR}/results/VCF
VCFOUT=${OUTDIR}/chr${CHR}.1000G.vcf.gz

# Filter 1000G
bcftools annotate -x INFO,FORMAT --rename-chrs ${SLURM_SUBMIT_DIR}/chr_names.txt $VCFIN |\
    bcftools norm -m +both |\
    bcftools view --genotype ^miss --types snps -m2 -M2 --min-af 0.01:minor |\
    bgzip -c > $VCFOUT

tabix -p vcf $VCFOUT
