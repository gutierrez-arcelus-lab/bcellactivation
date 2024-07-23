#!/usr/bin/bash

source /programs/biogrids.shrc
export PYTHON_X=2.7.2

INPUT=../data/rsid_23andMe.txt
OUTPUT=../data/rsid_23andMe_updated.txt

python /lab-share/IM-Gutierrez-e2/Public/References/dbSNP/liftRsNumber.py $INPUT > $OUTPUT
