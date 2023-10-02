source /programs/biogrids/biogrids.shrc
export PYTHON_X=2.7.2
export R_X=4.1

JUNC=./juncfiles.txt
OUT=./results_p1

# Intron clustering
python $HOME/leafcutter/clustering/leafcutter_cluster_regtools.py \
    -j $JUNC \
    -m 30 \
    -l 500000 \
    -p 0.00001 \
    -o ${OUT}/leafcutter

rm *.leafcutter.sorted.gz
