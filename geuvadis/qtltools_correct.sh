#!/bin/bash

source /programs/biogrids.shrc
export QTLTOOLS_X=1.3.1

BED=./ebv_phenotypes.bed
COV=./qtltools_cov.txt
OUT=./ebv_phenotypes_corrected.bed

QTLtools correct \
    --bed $BED \
    --cov $COV \
    --out $OUT \
    --normal

bgzip $OUT
tabix -p bed ${OUT}.gz
