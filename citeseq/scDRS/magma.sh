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

${TOOLS}/magma/magma --annotate window=50 --snp-loc $SNPLOC --gene-loc $GENELOC --out $OUT

# Magma gene analysis
REF=${TOOLS}/magma/aux_data/g1000_eur/g1000_eur
ANNOT=${OUT}.genes.annot
PVAL=./data/magma_input/${PREFIX}_snp_pval.tsv
OUTGENE=./output/${PREFIX}

${TOOLS}/magma/magma --bfile $REF --gene-annot $ANNOT --pval $PVAL N=10995 --out $OUTGENE

# scDRS munge
#ZSCORE=${OUTGENE}_zscore
#OUTGS=${OUTGENE}.gs
#
#scdrs munge-gs \
#    --out-file <out_file> \
#    --zscore-file <zscore_file> \
#    --weight zscore \
#    --n-max 1000


conda deactivate
