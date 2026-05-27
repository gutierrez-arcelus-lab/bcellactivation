#!/usr/bin/env bash

# ==============================================================================
# Description:  Compiles the statistical outputs from the differential splicing 
#               analysis into a structured `.Rdata` object. This object is 
#               required to interactively visualize the alternative splicing 
#               events using the LeafViz Shiny application.
# Input:        1. ./data/contrasts.txt (List of pairwise comparisons)
#               2. ./data/clusters_perind_numers.counts.gz (Cluster count matrix)
#               3. ./data/groups_${CONTRAST}.tsv (Sample-to-group mapping)
#               4. ./results/${CONTRAST}_cluster_significance.txt (P-values)
#               5. ./results/${CONTRAST}_effect_sizes.txt (dPSI values)
#               6. ./data/leafviz_annot (Genome/transcriptome annotation for UI)
# Output:       ./results/${CONTRAST}.Rdata (Ready for LeafViz upload)
# ==============================================================================

# ==============================================================================
# SLURM Configuration
# ==============================================================================

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=4gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute,bch-compute-pe
#SBATCH --array=1-6
#SBATCH --job-name=leafviz
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/leafviz.%A.%a.out
#SBATCH -e logs/leafviz.%A.%a.err

# ==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids (Specific to our institutional cluster).
# External users: Ensure R and the required LeafViz R packages are installed.
source /programs/biogrids/biogrids.shrc
export PYTHON_X=2.7.2
export R_X=4.1

cd "$SLURM_SUBMIT_DIR"

# ==============================================================================
# Input Parsing
# ==============================================================================

# Extract the specific contrast name for this array task (e.g., "BCRvs.unstday0")
CONTRAST=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' ./data/contrasts.txt )

# Define paths to the inputs required for this specific contrast
COUNT="./data/clusters_perind_numers.counts.gz"
GROUP="./data/groups_${CONTRAST}.tsv"

# ==============================================================================
# Compile LeafViz Object
# ==============================================================================

# Run the LeafViz preparation script.
# --code: A prefix/identifier used by the LeafViz UI to label this contrast.
# --meta_data_file: The 2-column condition mapping file used in the previous step.
# Positional arguments: 
#   1) Raw cluster counts
#   2) Cluster significance results
#   3) Effect size (dPSI) results
#   4) Directory containing the pre-built LeafViz annotation files
# -o: The final output .Rdata file.
# NOTE FOR EXTERNAL USERS: Update the script path below to point to your local 
# installation of prepare_results.R within the LeafCutter/leafviz directory.

/lab-share/IM-Gutierrez-e2/Public/tools/leafcutter/leafviz/prepare_results.R \
    --code "${CONTRAST}" \
    --meta_data_file "$GROUP" \
    "$COUNT" \
    "./results/${CONTRAST}_cluster_significance.txt" \
    "./results/${CONTRAST}_effect_sizes.txt" \
    "./data/leafviz_annot" \
    -o "./results/${CONTRAST}.Rdata"
