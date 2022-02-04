#!/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=8gb
#SBATCH --time=72:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=ADMIXTURE_refpanel
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/ADMIXTURE_refpanel.%j

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

VCF1000G=${SLURM_SUBMIT_DIR}/results/allchr.merged.pruned.1000G.vcf.gz
SAMPLES=${SLURM_SUBMIT_DIR}/refpanel_ids.txt
PREFIX=${SLURM_SUBMIT_DIR}/results/allchr.refpanel
VCFOUT=${PREFIX}.vcf.gz
BEDOUT=${PREFIX}.bed

# Extract data for ref panel samples from 1000G data
bcftools view --samples-file $SAMPLES --force-samples -O z -o $VCFOUT $VCF1000G

# Make plink files
plink --vcf $VCFOUT \
    --keep-allele-order \
    --make-bed \
    --out $PREFIX

# Run admixture
cd ${SLURM_SUBMIT_DIR}/results

K=3
admixture $BEDOUT $K -j4

# Run projection
MGB=allchr.merged.pruned.MGB
BEDMGB=${SLURM_SUBMIT_DIR}/results/${MGB}.bed
cp ${PREFIX}.${K}.P ${MGB}.${K}.P.in

admixture -P $BEDMGB $K -j4

