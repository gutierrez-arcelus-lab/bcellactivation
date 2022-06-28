#!/usr/bin/env bash

#ANNOT=/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.primary_assembly.annotation.gtf.gz
#EXON=./grch38_exon.txt.gz
ANNOT=/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz
EXON=./hg19_exon.txt.gz

$HOME/leafcutter/scripts/gtf_to_exons.R $ANNOT $EXON

