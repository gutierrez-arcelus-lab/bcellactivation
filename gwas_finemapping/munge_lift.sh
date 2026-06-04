#!/usr/bin/bash

source /programs/biogrids.shrc
export BCFTOOLS_X=1.22

LIFT_DATA="/lab-share/IM-Gutierrez-e2/Public/References/Genomes/hsapiens/lift_data"

GENOME37="${LIFT_DATA}/GRCh37.p13.genome.fa"
GENOME38="${LIFT_DATA}/GRCh38.p14.genome.fa"
CHAIN="${LIFT_DATA}/hg19ToHg38.over.chain.gz"
HEADER="./data/munge_col_header.tsv"

GWAS="langefeld"

SUMSTATS="./data/${GWAS}_munge.tsv"
VCF19="./data/${GWAS}_gwas_grch37.vcf.gz"
VCF38="./data/${GWAS}_gwas_grch38.vcf.gz"

bcftools +munge -C "$HEADER" -f "$GENOME37" -s stats -Ou "$SUMSTATS" |
    bcftools norm -c wx -f "$GENOME37" -Ou |
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Ou |
    bcftools sort -O z -o "$VCF19" -

bcftools index -t -f "$VCF19" 

# Lift over for plotting only with ATAC-seq tracks (Fig 5)
bcftools +liftover -Ou "$VCF19" -- -s "$GENOME37" -f "$GENOME38" -c "$CHAIN" --reject-type v --reject "./data/failed_to_lift_${GWAS}.vcf" |
    bcftools norm -f "$GENOME38" -Ou |
    bcftools sort -O z -o "$VCF38"

bcftools index -t -f "$VCF38" 
