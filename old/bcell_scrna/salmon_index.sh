#!/usr/bin/env bash

source /programs/biogrids.shrc
export SALMON_X=1.5.1

GENCODE=$HOME/lupus/data/gencode.transcripts.fa
EBV=$HOME/lupus/data/ebv_NC_007605.1_genome.fa
FASTA=$HOME/lupus/indices/salmon_human_ebvdecoy.fasta
DECOY=$HOME/lupus/indices/decoys.txt
OUT=$HOME/lupus/indices/SALMON_HUMAN_EBVDECOY

grep "^>" $EBV |\
    cut -d " " -f 1 |\
    sed -e 's/>//g' > $DECOY

cat $GENCODE $EBV > $FASTA

salmon index -t $FASTA -d $DECOY -p 8 -i $OUT 

rm $FASTA $DECOY
