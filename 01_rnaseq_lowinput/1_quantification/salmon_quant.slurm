#!/usr/bin/bash

# ==============================================================================
# Description:  Slurm array script to quantify transcript abundance using Salmon.
#               Handles multiple technical replicates/lanes per sample by parsing 
#               comma-separated FASTQ lists and passing them to Salmon.
# Input:        metadata.tsv (Headerless: ID, FQ1_list, FQ2_list)
#               Pre-built Salmon Index
# Output:       Transcript quantification directories for each sample.
# Dependencies: Biogrids (Salmon v1.10.1)
# ==============================================================================

# ------------------------------------------------------------------------------
# SLURM DIRECTIVES
# ------------------------------------------------------------------------------

#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=24gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-152
#SBATCH -J Salmon
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/Salmon-quant.%A-%a.out
#SBATCH -e logs/Salmon-quant.%A-%a.err

# ------------------------------------------------------------------------------
# ENVIRONMENT SETUP
# ------------------------------------------------------------------------------
# Note: Biogrids is our institutional software environment manager. 
# If running this script on a different HPC cluster, replace the following 
# two lines with your local equivalent (e.g., `conda activate salmon`).
source /programs/biogrids.shrc
export SALMON_X=1.10.1

# Change to directory where script was submitted from
cd ${SLURM_SUBMIT_DIR}

# ------------------------------------------------------------------------------
# INPUT PARSING & ARRAY BUILDING
# ------------------------------------------------------------------------------
META="./data/metadata.tsv"

# Extract Sample ID and comma-separated FASTQ names for the current array task
SAMPLEID=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $META)
FQ1=$(awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $META)
FQ2=$(awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $3 }' $META)

# Convert the comma-separated strings into Bash arrays by replacing commas with spaces
FQLIST1=(${FQ1//,/ })
FQLIST2=(${FQ2//,/ })

# Prepend the absolute path of the working directory to every element in the arrays
FQLIST1=( "${FQLIST1[@]/#/${TEMP_WORK}/fastq/lowinput/}" )
FQLIST2=( "${FQLIST2[@]/#/${TEMP_WORK}/fastq/lowinput/}" )

# Define Transcriptome Index path
INDEX="./index"

# ------------------------------------------------------------------------------
# OUTPUT DIRECTORY
# ------------------------------------------------------------------------------
OUT="./results/salmon/${SAMPLEID}"
mkdir -p "$OUT"

# ------------------------------------------------------------------------------
# EXECUTION
# ------------------------------------------------------------------------------
# Run Salmon quantification
# -l A            : Automatically infer the library type (e.g., strandedness)
# --seqBias       : Correct for sequence-specific biases at the 5' and 3' ends
# --posBias       : Correct for positional biases (e.g., 3' coverage bias in RNA-seq)
# --gcBias        : Correct for fragment-level GC content bias
# --numBootstraps : Run 20 bootstraps to estimate variance in abundance estimates 
salmon quant -i "$INDEX" -l A \
    -1 "${FQLIST1[@]}" -2 "${FQLIST2[@]}" \
    -p "$SLURM_CPUS_PER_TASK" \
    -o "$OUT" \
    --seqBias --posBias --gcBias \
    --numBootstraps 20
