#!/usr/bin/bash

# ==============================================================================
# Description:  Performs adapter trimming and quality control for paired-end 
#               RNA-seq data using Trim Galore and FastQC.
# Input:        A 3-column TSV named `metadata.tsv` in the `data/` directory.
#               Format: [sample_id, fastq_R1, fastq_R2] (long format, with header).
#               Each row represents a single pair of FASTQ files.
# Output:       Trimmed FASTQ files and FastQC reports in the output directory.
# ==============================================================================

# ==============================================================================
# SLURM Configuration
# ==============================================================================

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-115
#SBATCH --job-name=FASTQC
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/FASTQC.%A-%a.log

# ==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids (Specific to our institutional cluster).
# External users: Ensure Trim Galore, Cutadapt, and FastQC are installed and 
# available in your $PATH.
source /programs/biogrids.shrc
export FASTQC_X=0.11.9

# ==============================================================================
# Input Parsing
# ==============================================================================
# Define the metadata file path
META=${SLURM_SUBMIT_DIR}/data/metadata.tsv

# Parse the R1 and R2 fastq paths for the current array task.
# `awk 'FNR>1'` skips the header row. 
# `awk -v ARRID... 'FNR==ARRID'` matches the Slurm array ID to the data row.
FQ1=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }')
FQ2=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $3 }')

# ==============================================================================
# Output Definition
# ==============================================================================
# NOTE FOR EXTERNAL USERS: Update this path to your desired output directory.
OUT=${TEMP_WORK}/fastq/highinput
mkdir -p $OUT

# ==============================================================================
# Quality Control & Trimming
# ==============================================================================
# Run Trim Galore:
# --paired: Trims paired-end files synchronously
# --length 75: Discards any reads that become shorter than 75bp after trimming
# --fastqc: Automatically runs FastQC on the resulting trimmed files
trim_galore --paired \
    --length 75 \
    --fastqc \
    --cores $SLURM_CPUS_PER_TASK \
    -o $OUT \
    $FQ1 $FQ2
