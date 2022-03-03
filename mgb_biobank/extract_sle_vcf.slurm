#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=01:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-23
#SBATCH --job-name=SLEVCF
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /temp_work/ch229163/log/SLEVCF.%A.%a 

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

if [ "$SLURM_ARRAY_TASK_ID" -lt 23 ]; then

    CHR=chr${SLURM_ARRAY_TASK_ID}

elif [ "$SLURM_ARRAY_TASK_ID" -eq 23 ]; then

    CHR=chrX

fi

# convert Hg19 coordinations to Hg38
LABSHARE=/lab-share/IM-Gutierrez-e2/Public 
BED=${SLURM_SUBMIT_DIR}/sle_variants/sle_variants.bed 
CHAIN=/reference_databases/ReferenceGenome/liftover_chain/hg19/hg19ToHg38.over.chain.gz

PREFIX="${BED%.*}"
BED38=${PREFIX}_hg38.bed
FAIL=${PREFIX}_failtolift.bed

if [ ! -f "$BED38" ]; then
  
    liftOver $BED $CHAIN $BED38 $FAIL 
fi

POS=${PREFIX}_hg38.txt
awk -F'\t' '{ printf ("%s\t%s\n", $1, $2) }' $BED38 > $POS

# Extract SLE variants from MGB VCFs
for BATCH in {01..08} 10
do
    VCFIN=${LABSHARE}/mgb_biobank/04${BATCH}/vcf/${CHR}.dose.vcf.gz 
    VCFOUT=${TEMP_WORK}/VCF/${CHR}.04${BATCH}.sle.vcf.gz
    SAMPLES=./results/eur_females_04${BATCH}.txt

    bcftools view -R $POS --samples-file $SAMPLES --force-samples $VCFIN |\
    bcftools annotate -x INFO,FORMAT -O z -o $VCFOUT -

    tabix -p vcf $VCFOUT

    VCFSCHR+=($VCFOUT)
done

MERGEDCHR=${SLURM_SUBMIT_DIR}/results/VCF/${CHR}.sle.merged.vcf.gz

bcftools merge -O z -o $MERGEDCHR $( echo "${VCFSCHR[*]}" )
tabix -p vcf $MERGEDCHR
