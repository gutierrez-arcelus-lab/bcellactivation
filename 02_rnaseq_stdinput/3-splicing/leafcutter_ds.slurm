#!/usr/bin/env bash

# ==============================================================================
# Description:  Performs differential intron excision analysis between specific 
#               biological conditions. It utilizes a Dirichlet-multinomial model 
#               via the `leafcutter_ds.R` script to identify alternative splicing 
#               events (cluster shifts) driven by stimulation.
# Input:        1. ./data/contrasts.txt (List of pairwise comparisons)
#               2. ./data/clusters_perind_numers.counts.gz (Cluster count matrix)
#               3. ./data/groups_${CONTRAST}.tsv (Sample-to-group mapping)
#               4. ./data/exon.txt.gz (Genomic exon annotations for context)
# Output:       Differential splicing statistics and effect sizes per cluster 
#               in the ./results/${CONTRAST}/ directory.
# ==============================================================================

# ==============================================================================
# SLURM Configuration
# ==============================================================================

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-6
#SBATCH --job-name=leafcut
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/leafcutter_ds.%A.%a.out
#SBATCH -e logs/leafcutter_ds.%A.%a.err

# ==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids (Specific to our institutional cluster).
# External users: Ensure R and the required LeafCutter R packages are installed.

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
EXON="./data/exon.txt.gz"

# ==============================================================================
# Differential Splicing Analysis
# ==============================================================================
# Run the LeafCutter differential splicing R script.
# -e: Path to the exon definitions (adds gene/exon context to the output clusters).
# -p: Number of threads to use for parallel processing.
# -o: Prefix for the output results files.
# Positional arguments: 1) The count matrix, 2) The contrast group file.
# NOTE FOR EXTERNAL USERS: Update the script path below to point to your local 
# installation of leafcutter_ds.R.
/lab-share/IM-Gutierrez-e2/Public/tools/leafcutter/scripts/leafcutter_ds.R \
    -e "$EXON" \
    -p "$SLURM_CPUS_PER_TASK" \
    -o "./results/${CONTRAST}" \
    "$COUNT" \
    "$GROUP"
