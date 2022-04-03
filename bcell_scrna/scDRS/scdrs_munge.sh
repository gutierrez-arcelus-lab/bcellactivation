#!/usr/bin/bash

ZSCORE=./results/Bentham_zscore.tsv
OUT=./results/Bentham_zscore.gs

scdrs munge_gs $OUT \
    --zscore_file $ZSCORE \
    --weight zscore \
    --n_max 1000
