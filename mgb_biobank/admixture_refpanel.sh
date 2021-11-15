#!/bin/bash

DIR=./results
K=5

admixture $DIR/allchr.refpanel.bed $K -j4

mv ./allchr.refpanel.$K.P ./allchr.refpanel.$K.Q $DIR/
