#!/usr/bin/bash

GWAS=Bentham
GENES=./data/${GWAS}_gene_ids.txt
OUTDIR=./data/qtls_susie/${GWAS} 

mkdir -p ${OUTDIR}

awk 'FNR>1 { print $2 }' ./data/${GWAS}_genes.tsv > $GENES

for i in $( find /lab-share/IM-Gutierrez-e2/Public/External_datasets/eQTLcatalogue/susie -name "*lbf*" )
do
    zcat $i | grep -F -f $GENES | bgzip > ${OUTDIR}/$( basename $i ) 
done
