#!/usr/bin/env bash

source /programs/biogrids/biogrids.shrc
export R_X=4.1

ANNOT=/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v41.primary_assembly.annotation.gtf.gz
EXON=./data/exon.txt.gz

/lab-share/IM-Gutierrez-e2/Public/tools/leafcutter/scripts/gtf_to_exons.R $ANNOT $EXON
