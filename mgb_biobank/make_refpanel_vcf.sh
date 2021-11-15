#!/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.12

# Extract data for ref panel samples from 1000G data
OUTDIR=./results 
VCF1000G=$OUTDIR/allchr.1000G.vcf.gz
SAMPLES=./ref_panel_ids.txt
VCFOUT=$OUTDIR/allchr.refpanel.vcf.gz

bcftools view --samples-file $SAMPLES -O z -o $VCFOUT $VCF1000G

plink --vcf $VCFOUT \
    --keep-allele-order \
    --make-bed \
    --out $OUTDIR/allchr.refpanel

