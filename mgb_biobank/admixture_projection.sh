#!/bin/bash

DIR=/lab-share/IM-Gutierrez-e2/Public/vitor/ase/mgb_biobank
K=5

cd $DIR/results

cp allchr.refpanel.$K.P allchr.MGB.$K.P.in 

admixture -P -j4 allchr.MGB.bed $K | tee admix.MGBprojection.$K.log
