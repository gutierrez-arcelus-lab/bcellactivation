#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-23
#SBATCH --job-name=filter_1000G
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/filter_1000G.%A.%a

# Load tools
source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Chromosome
if [ "$SLURM_ARRAY_TASK_ID" -lt 23 ]; then

    CHR=$SLURM_ARRAY_TASK_ID

elif [ "$SLURM_ARRAY_TASK_ID" -eq 23 ]; then

    CHR=X

fi

# Define IO
LABSHARE=/lab-share/IM-Gutierrez-e2/Public
INDIR=${LABSHARE}/vcf_1000G
OUTDIR=${SLURM_SUBMIT_DIR}/VCF
VCFIN=${INDIR}/ALL.chr${CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
VCFOUT=${OUTDIR}/chr${CHR}.1000G.vcf.gz
SAMPLES=${SLURM_SUBMIT_DIR}/geuvadis_kgp_ids.txt

mkdir -p $OUTDIR

# Filter 1000G
bcftools view --samples-file $SAMPLES $VCFIN |\
    bcftools view --genotype ^miss --min-af 0.05:minor -m2 -M2 |\
    bcftools annotate -x INFO,FORMAT --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' - |\
    bgzip -c > $VCFOUT

tabix -p vcf $VCFOUT
