#!/usr/bin/bash

# ==============================================================================
# SLURM Configuration
# ==============================================================================

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=16gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-74
#SBATCH --job-name=fastqc
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/fastqc.%A-%a.log
#SBATCH -e logs/fastqc.%A-%a.err

# NOTE FOR EXTERNAL USERS: 
# 1. Change the array number below to match the number of samples (data rows) in metadata.tsv
# 2. Update the -o and -e log paths to a directory that exists on your system
#SBATCH --array=1-74
#SBATCH -o /temp_work/ch229163/log/fastqc.%A-%a.log
#SBATCH -e /temp_work/ch229163/log/fastqc.%A-%a.err

#==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids (Specific to our institutional cluster). 
# External users: Ensure Trim Galore, Cutadapt, FastQC, and Python are in your $PATH.
source /programs/biogrids.shrc
export CUTADAPT_X=4.1
export FASTQC_X=0.11.9
export PYTHON_X=3.8.8

# ==============================================================================
# 1. Input Parsing from Metadata
# ==============================================================================
# Shift array ID by 1 to account for the header row in the TSV
TARGET_ROW=$((SLURM_ARRAY_TASK_ID + 1))
META=${SLURM_SUBMIT_DIR}/data/metadata.tsv

# Extract sample information for the current array task
DONORID=$( awk -v row="$TARGET_ROW" 'NR==row { print $1 }' $META )
REPID=$( awk -v row="$TARGET_ROW" 'FNR==row { print $2 }' $META )
STIM=$( awk -v row="$TARGET_ROW" 'FNR==row { print $3 }' $META )
FQS1=$( awk -v row="$TARGET_ROW" 'FNR==row { print $4 }' $META )
FQS2=$( awk -v row="$TARGET_ROW" 'FNR==row { print $5 }' $META )

# ==============================================================================
# 2. Output Definition & FastQ Merging
# ==============================================================================
# NOTE FOR EXTERNAL USERS: Replace $TEMP_WORK with your desired output directory
OUTDIR=${TEMP_WORK}/fastq/highinput_merged
mkdir -p $OUTDIR

# Convert comma-separated string of files (from setup.R) into a bash array
FQLIST1=(${FQS1//,/ })
FQLIST2=(${FQS2//,/ })

# Define the names for the merged temporary fastq files
FQMERGE1=${OUTDIR}/${DONORID}_${REPID}_${STIM}_R1.fastq.gz
FQMERGE2=${OUTDIR}/${DONORID}_${REPID}_${STIM}_R2.fastq.gz

# Concatenate all fragments/resequenced lanes into single R1 and R2 files
cat ${FQLIST1[@]} > $FQMERGE1
cat ${FQLIST2[@]} > $FQMERGE2

# ==============================================================================
# 3. Quality Control & Trimming
# ==============================================================================
# Run Trim Galore (which internally runs Cutadapt and FastQC)
trim_galore --paired \
    --length 50 \
    --fastqc \
    --cores $SLURM_CPUS_PER_TASK \
    -o $OUTDIR \
    $FQMERGE1 $FQMERGE2

# Clean up the massive temporary merged files to save storage space
rm $FQMERGE1 $FQMERGE2
