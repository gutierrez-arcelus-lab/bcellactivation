#!/bin/bash

K=5

cd ./results

#cp allchr.1000G.$K.P allchr.MGB.$K.P.in 

admixture -P -j8 allchr.MGB.bed $K
