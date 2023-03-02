#!/usr/bin/bash

source /programs/biogrids.shrc

CITESEQ=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq/SN0268787/broad/hptmp/curtism/bwh10x/KW10456_Maria/221101_10X_KW10456_bcl/cellranger-6.1.1/GRCh38/BRI-1984/outs
BAM=${CITESEQ}/possorted_genome_bam.bam
BARCODES=${CITESEQ}/filtered_feature_bc_matrix/barcodes.tsv.gz
VCF=../../data/allchr.mgb.generegions.vcf.gz
OUT=./demuxlet_1984_results
LOG=${OUT}.log

demuxlet --sam $BAM \
    --vcf $VCF --field GP \
    --group-list $BARCODES \
    --min-BQ 20 \
    --min-uniq 10 \
    --alpha 0 --alpha 0.5 \
    --out $OUT 2> $LOG
