#!/usr/bin/bash

REF=./data/ref_1000G/g1000_eur
ANNOT=./results/Bentham.genes.annot
PVAL=./data/bentham_pval.tsv
OUT=./results/Bentham

$HOME/bin/magma/magma --bfile $REF --gene-annot $ANNOT --pval $PVAL N=10995 --out $OUT
