#!/usr/bin/bash

source /programs/biogrids.shrc

CITESEQ=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq/SN0273471/broad/hptmp/sgurajal/bwh10x/KW10598_mgutierrez/221211_10X_KW10598_bcl/cellranger-6.1.1/GRCh38/BRI-1988_hashing/outs
BAM=${CITESEQ}/possorted_genome_bam.bam
BARCODES=${CITESEQ}/filtered_feature_bc_matrix/barcodes.tsv.gz
VCF=./data/allchr.mgb.generegions.vcf.gz
OUT=./demuxlet_1988_results
LOG=${OUT}.log

demuxlet --sam $BAM \
    --vcf $VCF --field GP \
    --group-list $BARCODES \
    --min-BQ 20 \
    --min-uniq 10 \
    --alpha 0 --alpha 0.5 \
    --out $OUT 2> $LOG
