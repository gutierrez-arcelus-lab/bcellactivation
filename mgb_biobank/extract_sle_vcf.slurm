#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=01:00:00
#SBATCH -p bch-compute
#SBATCH --array=0-22
#SBATCH --job-name=SLEVCF
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/SLEVCF.%A.%a 

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

LABSHARE=/lab-share/IM-Gutierrez-e2/Public 

# Chromosome
ALLCHRS=( $(seq 1 22; echo X) )
CHR=chr${ALLCHRS[$SLURM_ARRAY_TASK_ID]}
SLE=${SLURM_SUBMIT_DIR}/sle_variants/sle_variants_hg38.txt

for BATCH in {0401..0408} 0410; do

    VCFIN=${LABSHARE}/mgb_biobank/${BATCH}/vcf/${CHR}.dose.vcf.gz
    VCFTMP=${TEMP_WORK}/VCF/${CHR}.${BATCH}.sle.tmp.vcf.gz
    SAMPLES=${SLURM_SUBMIT_DIR}/results/eur_females_${BATCH}.txt

    bcftools view -R $SLE -S $SAMPLES --force-samples $VCFIN |\
        bcftools norm -m +both |\
        bcftools annotate -x INFO,FORMAT -O z -o $VCFTMP

    tabix -p vcf $VCFTMP
    VCFSCHR+=($VCFTMP)
done

NVCF=$(echo "${#VCFSCHR[@]}")
OUTDIR=${TEMP_WORK}/VCF
VCFOUT=${OUTDIR}/${CHR}_intersect_sle

bcftools isec -p $VCFOUT -n=$NVCF -c both -O z -o $VCFOUT "${VCFSCHR[@]}"

MERGEDCHR=${SLURM_SUBMIT_DIR}/results/VCF/${CHR}.sle.merged.vcf.gz
VCFMAX=$(( $NVCF - 1 ))
VCFISEC=( $(eval echo ${VCFOUT}/000{0..$VCFMAX}.vcf.gz) )

bcftools merge "${VCFISEC[@]}" |\
    bcftools annotate -x INFO,FORMAT --set-id '%CHROM:%POS:%REF:%ALT' -O z -o $MERGEDCHR  

tabix -p vcf $MERGEDCHR
