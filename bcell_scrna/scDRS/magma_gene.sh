#!/usr/bin/bash

REF=./data/ref_1000G/g1000_eur
#ANNOT=./results/Bentham.genes.annot
ANNOT=./results/SleRiskGenes.genes.annot
PVAL=./data/summ_stats/bentham_pval.tsv
#OUT=./results/Bentham
OUT=./results/SleRiskGenes

$HOME/bin/magma/magma --bfile $REF --gene-annot $ANNOT --pval $PVAL N=10995 --out $OUT
