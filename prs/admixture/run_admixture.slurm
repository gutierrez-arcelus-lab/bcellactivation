#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -J ADMIXTURE
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/ADMIXTURE.%j.out
#SBATCH -e /temp_work/ch229163/log/ADMIXTURE.%j.err

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

## Extract data for ref panel samples from 1000G data
DIR=${TEMP_WORK}/vcf/prs 
VCF1000G=${DIR}/allchr.1000G.vcf.gz
SAMPLES=${SLURM_SUBMIT_DIR}/data/kgp_refpanel_ids.txt
PREFIX=allchr.refpanel
VCFREF=${DIR}/${PREFIX}.vcf.gz
BEDREF=${DIR}/${PREFIX}.bed

bcftools view --threads $SLURM_CPUS_PER_TASK \
    --samples-file $SAMPLES --force-samples \
    -O z -o $VCFREF \
    $VCF1000G

# Extract data for MGB samples
MERGED=${DIR}/allchr.merged.pruned.vcf.gz
SAMPLESMGB=${SLURM_SUBMIT_DIR}/data/mgb_ids.txt
PREFIXMGB=allchr.mgb
VCFMGB=${DIR}/${PREFIXMGB}.vcf.gz
BEDMGB=${DIR}/${PREFIXMGB}.bed

#bcftools view --threads $SLURM_CPUS_PER_TASK \
#    --samples-file $SAMPLESMGB \
#    -O z -o $VCFMGB \
#    $MERGED

# Make plink files
plink --vcf $VCFREF \
    --keep-allele-order \
    --make-bed \
    --out $PREFIX

#plink --vcf $VCFMGB \
#    --keep-allele-order \
#    --make-bed \
#    --out $PREFIXMGB

# Run admixture
cd ${SLURM_SUBMIT_DIR}

K=5
admixture $BEDREF $K -j8

# Run projection
cp ${PREFIX}.${K}.P ${PREFIXMGB}.${K}.P.in

admixture -P $BEDMGB $K -j8
