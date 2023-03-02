#!/usr/bin/bash

source /programs/biogrids/biogrids.shrc

CITESEQ=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/CITEseq_pilot_2/SN0257788/broad/hptmp/curtism/bwh10x/KW10170_Maria/220617_10X_KW10170_bcl/cellranger-6.1.1/GRCh38/BRI-1743_hashing/outs
BAM=${CITESEQ}/possorted_genome_bam.bam
BARCODES=${CITESEQ}/filtered_feature_bc_matrix/barcodes.tsv.gz
VCF=../../data/pilot2_genotypes.vcf.gz
OUT=./demuxlet/demuxlet_pilot2
LOG=${OUT}.log

demuxlet --sam $BAM \
    --vcf $VCF --field GT \
    --group-list $BARCODES \
    --min-BQ 20 \
    --min-uniq 10 \
    --alpha 0 --alpha 0.5 \
    --out $OUT 2> $LOG
