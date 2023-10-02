#!/usr/bin/env bash

source /programs/biogrids.shrc
export GATK_X=4.1.4.1

BAM=/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/CITEseq_pilot_2/SN0257788/broad/hptmp/curtism/bwh10x/KW10170_Maria/220617_10X_KW10170_bcl/cellranger-6.1.1/GRCh38/BRI-1743_hashing/outs/possorted_genome_bam.bam
GENOME=/lab-share/IM-Gutierrez-e2/Public/References/Genomes/hsapiens/GRCh38/GRCh38.primary_assembly.genome.fa
VCF=./demuxlet/samples_allsnps_asereadc.vcf.gz
OUT=./demuxlet/asereadcounter

gatk ASEReadCounter \
    -R $GENOME \
    -I $BAM \
    -V $VCF \
    -O $OUT \
    --min-base-quality 20 \
    --min-mapping-quality 20
