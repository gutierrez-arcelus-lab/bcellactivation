#!/usr/bin/bash

# ==============================================================================
# Description:  Slurm array script for performing quality control and adapter 
#               trimming on paired-end RNA-seq data using Trim Galore & FastQC.
# Input:        description_hsapiens_longformat.tsv (List of samples)
#               - Column 8: Path to Read 1 (FASTQ)
#               - Column 9: Path to Read 2 (FASTQ)
# Output:       Trimmed FASTQ files and FastQC HTML reports.
# Dependencies: Biogrids (Trim Galore, FastQC v0.11.9)
# ==============================================================================

# ------------------------------------------------------------------------------
# SLURM DIRECTIVES
# ------------------------------------------------------------------------------

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-732
#SBATCH --job-name=FASTQC
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o /logs/FASTQC.%A-%a.out
#SBATCH -e /logs/FASTQC.%A-%a.err

# ------------------------------------------------------------------------------
# ENVIRONMENT SETUP
# ------------------------------------------------------------------------------
# Note: Biogrids is our institutional software environment manager. 
# If running this script on a different HPC cluster, replace the following 
# two lines with your local equivalent (e.g., `module load fastqc/0.11.9` 
# or `conda activate rna_seq_env`).
source /programs/biogrids.shrc
export FASTQC_X=0.11.9

# Change to the directory where the script was submitted from
cd ${SLURM_SUBMIT_DIR}

# ------------------------------------------------------------------------------
# INPUT PARSING
# ------------------------------------------------------------------------------
# Define the sample sheet path
SAMPLELIST="./data/metadata_longformat.tsv"

# Extract the input FASTQ paths for the current array task.
# FNR>1 skips the header row.
# FNR==ARRID matches the Slurm Array ID to the correct row in the TSV.
# Read 1 is extracted from column 2, Read 2 is extracted from column 3.
FQ1=$( awk 'FNR>1' $SAMPLELIST | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }')
FQ2=$( awk 'FNR>1' $SAMPLELIST | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $3 }')

# ------------------------------------------------------------------------------
# OUTPUT DIRECTORY
# ------------------------------------------------------------------------------

# Adjust to your desired output directory
OUT="./results"
mkdir -p $OUT

# ------------------------------------------------------------------------------
# EXECUTION
# ------------------------------------------------------------------------------
# Run Trim Galore with FastQC integration
# --length 50 : Discard reads that become shorter than 50bp after trimming
# --fastqc    : Run FastQC in the default mode on the FASTQ file once trimming is complete
trim_galore --paired \
    --length 50 \
    --fastqc \
    --cores $SLURM_CPUS_PER_TASK \
    -o $OUT \
    $FQ1 $FQ2
