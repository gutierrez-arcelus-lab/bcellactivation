#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=ADMIXTURE_refpanel
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /home/ch229163/log/ADMIXTURE_refpanel.%j

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

VCF1000G=${SLURM_SUBMIT_DIR}/results/allchr.merged.pruned.1000G.vcf.gz
SAMPLES=${SLURM_SUBMIT_DIR}/refpanel_ids.txt
PREFIX=allchr.refpanel 
OUT=${SLURM_SUBMIT_DIR}/results/${PREFIX} 
VCFOUT=${OUT}.vcf.gz
BEDOUT=${OUT}.bed

# Extract data for ref panel samples from 1000G data
bcftools view --samples-file $SAMPLES --force-samples -O z -o $VCFOUT $VCF1000G

# Make plink files
plink --vcf $VCFOUT \
    --keep-allele-order \
    --make-bed \
    --out $OUT

# Run admixture
K=4
admixture $BEDOUT $K -j4

mv ${SLURM_SUBMIT_DIR}/${PREFIX}.${K}.P ${SLURM_SUBMIT_DIR}/${PREFIX}.${K}.Q ${SLURM_SUBMIT_DIR}/results/
