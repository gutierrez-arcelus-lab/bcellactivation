#!/usr/bin/bash

H5=./data/expression/bcells_singlet_seurat.h5ad
#GS=./results/Bentham_zscore.gs
#OUT=./results
GS=./results/Bentham_zscore_scDRSpaper.gs
OUT=./results/Bentham_data_scDRSpaper

mkdir -p $OUT

scdrs compute_score $H5 human $GS human $OUT \
    --flag_filter_data False \
    --flag_raw_count False \
    --n_ctrl 1000 \
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score True
