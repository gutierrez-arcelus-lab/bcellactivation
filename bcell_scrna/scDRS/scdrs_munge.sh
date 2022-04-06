#!/usr/bin/bash

#ZSCORE=./results/Bentham_zscore.tsv
#OUT=./results/Bentham_zscore.gs

ZSCORE=./data/Bentham_zscore_from_scDRS.tsv
OUT=./results/Bentham_zscore_scDRSpaper.gs

scdrs munge_gs $OUT \
    --zscore_file $ZSCORE \
    --weight zscore \
    --n_max 1000
