#!/usr/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.22
export PLINK2_X=2.00a7

LIFT_DATA="/lab-share/IM-Gutierrez-e2/Public/References/Genomes/hsapiens/lift_data"
GENOME37="${LIFT_DATA}/GRCh37.p13.genome.fa"

VCF="/reference_databases/1000G_VCF/phase3/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
IDS="./data/tgp_eur.txt"
OUT="./data/chr1_TNFSF4"

# VCF processing
bcftools view -r 1:172191475-174191475 --samples-file "$IDS" --force-samples -Ou "$VCF" |
    bcftools norm -m -any -Ou |
    bcftools view -i 'MAC>=2' -Ou |
    bcftools annotate --rename-chrs <(echo "1 chr1") -Ou |
    bcftools norm -c wx -f "$GENOME37" -Ou | 
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Ou |
    bcftools sort -O z -o "${OUT}.vcf.gz"

# LD
plink2 --vcf "${OUT}.vcf.gz" --r-unphased square ref-based --memory 8000 --out "${OUT}_r2"
