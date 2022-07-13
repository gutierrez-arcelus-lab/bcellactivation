#!/usr/bin/env bash

ANNOT38=/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.primary_assembly.annotation.gtf.gz
EXON38=./grch38_exon.txt.gz
ANNOT19=/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz
EXON19=./hg19_exon.txt.gz

$HOME/leafcutter/scripts/gtf_to_exons.R $ANNOT38 $EXON38
$HOME/leafcutter/scripts/gtf_to_exons.R $ANNOT19 $EXON19

