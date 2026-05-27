#!/usr/bin/bash

# ==============================================================================
# Description:  Re-quantifies transcript abundance using the pooled metadata sheet. 
#               Crucially, this script maps against the GENCODE v38 index 
#               generated during the low-input pipeline to allow for a direct, 
#               fold-change comparison between the two datasets 
# Input:        ./data/metadata.tsv (Pooled sample sheet)
#               Trimmed FASTQ files (paths provided via metadata)
# Output:       Transcript quantification files (`quant.sf`) and bootstraps 
#               in sample-specific subdirectories inside `./results/salmon/`.
# ==============================================================================

# ==============================================================================
# SLURM Configuration
# ==============================================================================

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=8gb
#SBATCH --time=02:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-62
#SBATCH --job-name=Salmon
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/salmon.%A-%a.out
#SBATCH -e logs/salmon.%A-%a.err

# ==============================================================================
# Environment Setup
# ==============================================================================

# Load BioGrids (Specific to our institutional cluster).
# External users: Ensure Salmon is installed and available in your $PATH.
source /programs/biogrids.shrc
export SALMON_X=1.10.1

# ==============================================================================
# Input Parsing
# ==============================================================================
META="${SLURM_SUBMIT_DIR}/data/metadata.tsv"

# Parse the sample ID and comma-separated FASTQ lists.
# `awk 'FNR>1'` skips the header row.
SAMPLEID=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $1 }' )
FQ1=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' ) 
FQ2=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $3 }' ) 

# Convert the comma-separated strings into bash arrays by replacing commas with spaces.
# This allows Salmon to merge the technical replicates pooled in the previous step.
FQLIST1=( $(echo "$FQ1" | tr "," " ") )
FQLIST2=( $(echo "$FQ2" | tr "," " ") )

# ==============================================================================
# Cross-Pipeline Index Reference
# ==============================================================================
# Directs Salmon to the GENCODE v38 index built in the low-input pipeline, 
# ensuring strict annotation version parity for downstream fold-change comparisons.
INDEX="../../01_rnaseq_lowinput/1_quantification/index"

# ==============================================================================
# Output Definition
# ==============================================================================
OUT="${SLURM_SUBMIT_DIR}/results/salmon/${SAMPLEID}"
mkdir -p "$OUT"

# ==============================================================================
# Salmon Quantification
# ==============================================================================
# Run Salmon:
# -i: The shared decoy-aware v38 gentrome index
# -l A: Automatically infer the library type (e.g., strandedness)
# -1 / -2: Pass the bash arrays of R1 and R2 fastq files
# Bias flags: Corrects for sequence, positional, and GC biases.
# --numBootstraps 20: Computes 20 bootstrap samples to estimate technical 
#                     variance in abundance estimates.

salmon quant -i "$INDEX" -l A \
    -1 "${FQLIST1[@]}" -2 "${FQLIST2[@]}" \
    -p "$SLURM_CPUS_PER_TASK" \
    -o "$OUT" \
    --seqBias --posBias --gcBias \
    --numBootstraps 20
