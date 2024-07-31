#!/usr/bin/env bash

source /programs/biogrids.shrc
source ${HOME}/miniconda3/etc/profile.d/conda.sh

conda activate scDRS

# GWAS prefix
PREFIX=sle_bentham

# Magma annotate
TOOLS=/lab-share/IM-Gutierrez-e2/Public/tools
SNPLOC=./data/magma_input/${PREFIX}_snp_loc.tsv
GENELOC=./data/magma_input/gene_loc_hg19.tsv
OUT=./output/${PREFIX}_annot

${TOOLS}/magma/1.10/magma --annotate \
    window=50 \
    --snp-loc $SNPLOC \
    --gene-loc $GENELOC \
    --out $OUT

# Magma gene analysis
REF=${TOOLS}/magma/aux_data/g1000_eur/g1000_eur
ANNOT=${OUT}.genes.annot
PVAL=./data/magma_input/${PREFIX}_snp_pval.tsv
OUTGENE=./output/${PREFIX}

${TOOLS}/magma/1.10/magma \
    --bfile $REF \
    --gene-annot $ANNOT \
    --pval $PVAL \
    N=10995 \
    --out $OUTGENE

# scDRS munge
ZSCORE=${OUTGENE}_zscore.tsv
OUTGS=${OUTGENE}.gs

awk -v OFS='\t' '{print $1,$8}' ${OUTGENE}.genes.out | sed "1,1s/ZSTAT/$PREFIX/" > $ZSCORE

scdrs munge-gs \
    --out-file $OUTGS \
    --zscore-file $ZSCORE \
    --weight zscore \
    --n-max 1000

# scDRS compute_score
H5=./data/scdrs_input/bcells.h5ad
OUTSCDRS=./output

scdrs compute_score \
    --h5ad-file $H5 \
    --h5ad-species human \
    --gs-file $OUTGS \
    --gs-species human \
    --n-ctrl 1000 \
    --out-folder $OUTSCDRS \
    --flag_filter_data False \
    --flag_raw_count False \
    --flag_return_ctrl_raw_score False \
    --flag_return_ctrl_norm_score False

conda deactivate
