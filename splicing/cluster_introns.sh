#!/usr/bin/env bash

export PYTHON_X=2.7.2

echo -n '' > juncfiles.txt

for i in $(find ./mapping -name *.junc)
do
    echo $i >> juncfiles.txt
done

python $HOME/leafcutter/clustering/leafcutter_cluster_regtools.py \
    -j juncfiles.txt \
    -m 50 \
    -l 500000 \
    -o scharer 

rm *.junc.scharer.sorted.gz
