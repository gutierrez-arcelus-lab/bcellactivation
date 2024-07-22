#!/usr/bin/bash

source /programs/biogrids.shrc
source $HOME/miniconda3/etc/profile.d/conda.sh

conda activate vireo

DIR=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_citeseq/SN0268787/broad/hptmp/curtism/bwh10x/KW10456_Maria/221101_10X_KW10456_bcl/cellranger-6.1.1/GRCh38/BRI-1984/outs 
BAM=${DIR}/possorted_genome_bam.bam
BARCODE=${DIR}/filtered_feature_bc_matrix/barcodes.tsv.gz
VARIANTS=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/citeseq/demultiplexing/genotypes/data/allchroms_mgb_variants.vcf.gz
OUT=./data/1984

mkdir -p $OUT

cellsnp-lite -s $BAM -b $BARCODE -R $VARIANTS -O $OUT -p 10 --minCOUNT 20 --gzip --genotype

VCF=/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/citeseq/demultiplexing/genotypes/data/allchroms_mgb.vcf.gz
RES=./results/1984

mkdir -p $RES

vireo -c $OUT -d $VCF -o $RES --genoTag=GT

conda deactivate
