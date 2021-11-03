#!/bin/bash

source /programs/biogrids.shrc

DIR=/lab-share/IM-Gutierrez-e2/Public/vitor/mgb_biobank 
VCF=$DIR/allchr.vcf.gz
OUT=$DIR/mgb
LOG=$DIR/mgb.pca.log

QTLtools pca --vcf $VCF --center --scale --out $OUT --log $LOG
