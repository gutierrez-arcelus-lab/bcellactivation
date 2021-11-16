#!/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

DIR=/lab-share/IM-Gutierrez-e2/Public/vitor/ase/mgb_biobank/ 
VCF1000G=$DIR/results/allchr.1000G.vcf.gz
SAMPLES=$DIR/refpanel_ids.txt
PREFIX=allchr.refpanel 
OUT=$DIR/results/$PREFIX 
VCFOUT=$OUT.vcf.gz
BEDOUT=$OUT.bed

# Extract data for ref panel samples from 1000G data
bcftools view --samples-file $SAMPLES --force-samples -O z -o $VCFOUT $VCF1000G

# Make plink files
plink --vcf $VCFOUT \
    --keep-allele-order \
    --make-bed \
    --out $OUT

# Run admixture
K=3
admixture $BEDOUT $K -j4

mv $DIR/$PREFIX.$K.P $DIR/$PREFIX.$K.Q $DIR/results/
