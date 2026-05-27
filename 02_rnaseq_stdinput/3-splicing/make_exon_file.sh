#!/usr/bin/env bash

# ==============================================================================
# Description:  Parses the GENCODE GTF annotation file to generate a compressed 
#               exon definitions file. This file is used by LeafCutter during 
#               the differential splicing step to map anonymous splice junctions 
#               back to known genes and exons, providing biological context to 
#               the output clusters.
# Input:        ../1-mapping/data/gencode.v41.primary_assembly.annotation.gtf
# Output:       ./data/exon.txt.gz
# ==============================================================================

# ==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids (Specific to our intitutional HPC) 
# Note: LeafCutter's gtf_to_exons.R requires R to execute.
source /programs/biogrids/biogrids.shrc
export R_X=4.1

# ==============================================================================
# Variable Definition
# ==============================================================================

# Define the path to the exact same GTF used during the STAR mapping step
ANNOT="../1-mapping/data/gencode.v41.primary_assembly.annotation.gtf"

# Define the output file path
EXON="./data/exon.txt.gz"

# ==============================================================================
# Generate Exon Definitions
# ==============================================================================
# Run the LeafCutter GTF parsing script.
# NOTE FOR EXTERNAL USERS: Update the script path below to point to your local 
# installation of gtf_to_exons.R within the LeafCutter directory.

/lab-share/IM-Gutierrez-e2/Public/tools/leafcutter/scripts/gtf_to_exons.R "$ANNOT" "$EXON"
