#!/bin/bash

source /programs/biogrids.shrc
export QTLTOOLS_X=1.3.1

BED=./ebv_phenotypes.bed
#COV=./qtltools_cov.txt
COV=./qtltools_cov_v2.txt
#OUT=./ebv_phenotypes_corrected.bed
OUT=./ebv_phenotypes_corrected_v2.bed

QTLtools correct \
    --bed $BED \
    --cov $COV \
    --out $OUT \
    --normal

bgzip $OUT
tabix -p bed ${OUT}.gz
