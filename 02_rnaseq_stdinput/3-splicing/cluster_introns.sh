#!/usr/bin/bash

# ==============================================================================
# Description:  Performs intron clustering across all samples using LeafCutter.
#               This step aggregates overlapping split reads (junctions) into 
#               distinct clusters, which represent alternative splicing events 
#               to be statistically tested downstream.
# Input:        ./data/junction_files.txt (List of .junc files from setup.R)
# Output:       1. ./data/clusters_perind_numers.counts.gz (Count matrix)
#               2. ./data/clusters_perind.counts.gz (Used for diff splicing)
# ==============================================================================

# ==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids. This is our intitutional software manager. 
# Note: LeafCutter's legacy clustering script strictly relies on Python 2.7.
source /programs/biogrids/biogrids.shrc
export PYTHON_X=2.7.2
export R_X=4.1

# ==============================================================================
# Variable Definition
# ==============================================================================
# The text file containing all the paths to the individual .junc files
JUNC="./data/junction_files.txt"

# ==============================================================================
# Intron Clustering
# ==============================================================================
# Run LeafCutter clustering:
# -j: Text file containing paths to junction files.
# -p 0.01: Minimum fraction of reads in a cluster that must support a junction 
#          for it to be included (1% threshold removes rare, noisy splice junctions).
# -o: Prefix for the output count matrices.
# NOTE FOR EXTERNAL USERS: Ensure the python script path points to your local 
# LeafCutter installation directory.
python /lab-share/IM-Gutierrez-e2/Public/tools/leafcutter/clustering/leafcutter_cluster_regtools.py \
    -j $JUNC \
    -p 0.01 \
    -o "./data/clusters"

# ==============================================================================
# Cleanup
# ==============================================================================
# Remove intermediate sample-specific sorted cluster files to save disk space,
# as the downstream differential splicing steps only require the compiled matrices.
rm *.clusters.sorted.gz
