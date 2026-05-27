#!/usr/bin/bash

# ==============================================================================
# Description:  Quantifies transcript-level abundance using Salmon. 
# Input:        `metadata.tsv` (Format: sample_id, fastq_R1, fastq_R2). 
#               Note: Multiple fastq files for a single sample must 
#               be a comma-separated string.
# Output:       Transcript quantification files (`quant.sf`) in sample-specific 
#               subdirectories inside `./results/`.
# ==============================================================================

# ==============================================================================
# SLURM Configuration
# ==============================================================================

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-74
#SBATCH --job-name=Salmon
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/Salmon-quant.%A-%a.out
#SBATCH -e logs/Salmon-quant.%A-%a.err

# ==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids (Specific to our institutional cluster).
# External users: Ensure Salmon is installed and available in your $PATH.
source /programs/biogrids.shrc
export SALMON_X=1.10.1

cd "$SLURM_SUBMIT_DIR"

# ==============================================================================
# Input Parsing
# ==============================================================================
# Reference the metadata sheet used during the mapping step
META="../1-mapping/data/metadata.tsv"

# Parse the sample ID and comma-separated FASTQ lists
SAMPLEID=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $1 }' )
FQ1=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' ) 
FQ2=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $3 }' ) 

# Convert the comma-separated strings into bash arrays by replacing commas with spaces.
# This allows Salmon to merge multiple fastq files on the fly.
FQLIST1=( $(echo $FQ1 | tr "," " ") )
FQLIST2=( $(echo $FQ2 | tr "," " ") )

# Index
INDEX="./index"

# ==============================================================================
# Output Definition
# ==============================================================================
OUT="./results/${SAMPLEID}"
mkdir -p "$OUT"

# ==============================================================================
# Salmon Quantification
# ==============================================================================
# Run Salmon:
# -i: The decoy-aware gentrome index
# -l A: Automatically infer the library type (e.g., strandedness)
# -1 / -2: Pass the bash arrays of R1 and R2 fastq files
# Bias flags: Corrects for sequence-specific priming (--seqBias), 
#             read start position biases (--posBias), and GC content (--gcBias).

salmon quant -i "$INDEX" -l A \
    -1 "${FQLIST1[@]}" -2 "${FQLIST2[@]}" \
    -p "$SLURM_CPUS_PER_TASK" \
    -o "$OUT" \
    --seqBias --posBias --gcBias
