#! /usr/bin/env bash

source /programs/biogrids.shrc
export SAMTOOLS_X=1.13

LABSHR=/lab-share/IM-Gutierrez-e2/Public
SAMPLES=( $(cat ${LABSHR}/vitor/ase/bcell_bulk/bcell_samples.txt) )
BAMS=( "${SAMPLES[@]/#/$TEMP_WORK/gatk/}" )
BAMS=( "${BAMS[@]/%/_recal.bam}" )
GENOME=${LABSHE}/vitor/ase/data/GRCh38.primary_assembly.genome.fa
OUTDIR=${TEMP_WORK}/gatk
PREFIX=${OUTDIR}/MG8989
MERGED=${PREFIX}.merged.bam


samtools merge -@ 4 $MERGED "${BAMS[@]}"

