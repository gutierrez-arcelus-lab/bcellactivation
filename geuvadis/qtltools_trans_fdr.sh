#! /usr/bin/env bash

source /programs/biogrids.shrc
export R_X=4.1

OUTDIR=${SLURM_SUBMIT_DIR}/results/qtltools 
NOMOUT=${OUTDIR}/trans.nominal
SEED=123
OUTPERMUT=${OUTDIR}/trans.perm.${SEED}
OUT=${OUTDIR}/trans.fdr.txt


Rscript ${SLURM_SUBMIT_DIR}/qtltools_runFDR_ftrans.R \
    ${NOMOUT}.hits.txt.gz \
    ${OUTPERMUT}.hits.txt.gz \
    $OUT
