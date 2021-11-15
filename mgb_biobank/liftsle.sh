#!/bin/bash

BED=./sle_variants/sle.bed 
CHAIN=/reference_databases/ReferenceGenome/liftover_chain/hg19/hg19ToHg38.over.chain.gz
OUT=./sle_variants/sle_hg38.bed
FAIL=./sle_variants/sle_failtolift.bed

liftOver $BED $CHAIN $OUT $FAIL 
