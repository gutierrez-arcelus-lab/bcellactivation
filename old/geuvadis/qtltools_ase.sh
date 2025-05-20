#!/usr/bin/env bash

source /programs/biogrids.shrc

SAMPLE=ERR188022
LABSHR=/lab-share/IM-Gutierrez-e2/Public
BAM=${LABSHR}/vitor/ase/geuvadis/results/star/pass_wasp/${SAMPLE}.uniq.bam 
SORT=${LABSHR}/vitor/ase/geuvadis/results/star/pass_wasp/${SAMPLE}.uniq.sort.bam
VCF=${LABSHR}/vitor/ase/geuvadis/results/gatk/${SAMPLE}.knownSNVs.biallelic.het.vcf 
GENOME=${LABSHR}/vitor/ase/data/GRCh38.primary_assembly.genome.fixnames.fa 
GTF=${LABSHR}/vitor/ase/data/gencode.v38.primary_assembly.annotation.gtf
OUT=${LABSHR}/vitor/ase/geuvadis/results/qtltools/${SAMPLE}

#samtools sort $BAM > $SORT
#samtools index $SORT

QTLtools ase \
    --bam $SORT \
    --vcf $VCF \
    --ind $SAMPLE \
    --mapq 255  \
    --fasta $GENOME \
    --gtf $GTF \
    --cov 8 \
    --out $OUT
