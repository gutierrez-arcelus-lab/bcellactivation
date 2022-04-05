#!/usr/bin/bash

SNP_LOC=./data/summ_stats/bentham_snp_loc.tsv
GENE_LOC=./data/gencode_genes/genes.v19.hg19.loc
#GENE_LOC=./data/ncbi_genes/NCBI37.3.gene.loc
OUT=./results/Bentham

$HOME/bin/magma/magma --annotate --snp-loc $SNP_LOC --gene-loc $GENE_LOC --out $OUT 
