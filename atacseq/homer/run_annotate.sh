#!/usr/bin/bash

source /programs/biogrids.shrc
export HOMER_X=4.11

STIM=DN2

GENOME=/lab-share/IM-Gutierrez-e2/Public/References/Genomes/hsapiens/GRCh38/GRCh38.primary_assembly.genome.fa
BED=./data/${STIM}.bed
MOTIF=$(echo ./results/${STIM}/knownResults/known*.motif)
OUT=./results/${STIM}/annotate.txt

annotatePeaks.pl $BED $GENOME -mask -m ${MOTIF[@]} > $OUT
