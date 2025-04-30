#!/usr/bin/env bash

source /programs/biogrids/biogrids.shrc
export R_X=4.1

/lab-share/IM-Gutierrez-e2/Public/tools/leafcutter/leafviz/gtf2leafcutter.pl \
    -o data/leafviz_annot \
    /lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v41.primary_assembly.annotation.gtf
