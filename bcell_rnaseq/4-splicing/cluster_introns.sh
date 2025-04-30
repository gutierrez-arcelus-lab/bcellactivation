#!/usr/bin/bash

source /programs/biogrids/biogrids.shrc
export PYTHON_X=2.7.2
export R_X=4.1

JUNC=./data/junction_files.txt

# Intron clustering
python $HOME/leafcutter/clustering/leafcutter_cluster_regtools.py \
    -j $JUNC \
    -o ./data/clusters

rm *.clusters.sorted.gz
