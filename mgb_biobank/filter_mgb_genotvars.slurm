#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=198
#SBATCH --job-name=vcfgt
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/vcf_extract_genot.%A.%a

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Define directory and chromosome given array design file
DESIGN=${SLURM_SUBMIT_DIR}/array_design.txt
BATCH=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $DESIGN )
CHR=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $DESIGN )

# Define IO
LABSHARE=/lab-share/IM-Gutierrez-e2/Public
VCFIN=${LABSHARE}/mgb_biobank/$BATCH/vcf/chr${CHR}.dose.vcf.gz 
SAMPLES=${SLURM_SUBMIT_DIR}/results/females_${BATCH}.txt

OUTDIR=${TEMP_WORK}/VCF/genotyped
VCFOUT=${OUTDIR}/chr${CHR}.MGB.${BATCH}.vcf.gz

mkdir -p $OUTDIR

bcftools view --threads $SLURM_CPUS_PER_TASK --include '%ID=@/temp_work/ch229163/allchr_typed.txt' $VCFIN |\
    bcftools view --threads $SLURM_CPUS_PER_TASK --samples-file $SAMPLES --force-samples |\
    bcftools annotate --threads $SLURM_CPUS_PER_TASK -x INFO,FORMAT - |\
    bgzip -c > $VCFOUT    

tabix -p vcf $VCFOUT
