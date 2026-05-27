#!/usr/bin/env bash

# ==============================================================================
# Description:  Converts a standard GTF annotation file into the optimized 
#               format required by the LeafViz interactive Shiny app. 
# Input:        ../1-mapping/data/gencode.v41.primary_assembly.annotation.gtf
# Output:       Optimized annotation files in the `data/leafviz_annot/` directory.
# ==============================================================================

# ==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids (This is specific to our institutional cluster)
# Note: This specific LeafCutter utility is a Perl script, but relies on R 
# environment variables for downstream compatibility.

source /programs/biogrids/biogrids.shrc
export R_X=4.1

# ==============================================================================
# Generate LeafViz Annotation Database
# ==============================================================================
# Run the LeafCutter GTF parsing script.
# -o: The output directory where the optimized annotation files will be saved.
# Positional argument: The path to the exact same GTF used for STAR mapping.
# NOTE FOR EXTERNAL USERS: Update the script path below to point to your local 
# installation of gtf2leafcutter.pl within the LeafCutter directory.
#
/lab-share/IM-Gutierrez-e2/Public/tools/leafcutter/leafviz/gtf2leafcutter.pl \
    -o "data/leafviz_annot" \
    "../1-mapping/data/gencode.v41.primary_assembly.annotation.gtf"
