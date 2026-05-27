#!/usr/bin/bash

# ==============================================================================
# Description:  Builds a decoy-aware quasi-mapping index for Salmon. 
#               By providing the entire genome as a "decoy" background alongside 
#               the target transcriptome, Salmon can accurately discard reads 
#               that map better to unannotated/intergenic genomic regions rather 
#               than falsely assigning them to similar transcripts.
# Input:        1. GRCh38.primary_assembly.genome.fa.gz (Genome FASTA)
#               2. gencode.v41.primary_assembly.transcripts.fa (Transcript FASTA)
# Output:       Compiled Salmon index directory (`./index/`)
# ==============================================================================

# ==============================================================================
# SLURM Configuration
# ==============================================================================

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=24gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -J Salmon-index
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/Salmon-index.%j.out
#SBATCH -e logs/Salmon-index.%j.err

# ==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids (Specific to our institutional cluster).
# External users: Ensure Salmon is installed and available in your $PATH.
source /programs/biogrids.shrc
export SALMON_X=1.10.1

cd $SLURM_SUBMIT_DIR

# ==============================================================================
# Variable Definition
# ==============================================================================
GENOME="./data/GRCh38.primary_assembly.genome.fa.gz"
GENCODE="./data/gencode.v41.primary_assembly.transcripts.fa"

# Intermediate decoy/gentrome paths
GENTROME="./data/salmon/gentrome.fa.gz" 
DECOYS="./data/salmon/decoys.txt"

# Final index output directory
OUT="./index"

# Create required directories if they do not exist
mkdir -p ./data/salmon
mkdir -p "$OUT"

# ==============================================================================
# 1. Generate Decoys List
# ==============================================================================
# Extract the sequence headers (chromosome names) from the genome FASTA file.
# These identifiers tell Salmon which sequences in the gentrome are the background
# decoys rather than target transcripts.
grep "^>" <( gunzip -c $GENOME ) | cut -d " " -f 1 > $DECOYS

# Clean up the headers by removing the '>' character
sed -i.bak -e 's/>//g' $DECOYS

# ==============================================================================
# 2. Build the Gentrome
# ==============================================================================
# Concatenate the annotated transcripts and the whole genome into a single file.
cat $GENCODE $GENOME > $GENTROME

# ==============================================================================
# 3. Compile Salmon Index
# ==============================================================================
# Run Salmon indexing:
# -t: Target gentrome FASTA
# -i: Output directory for the compiled index
# --decoys: Text file listing the names of the decoy (genomic) sequences
# -k 31: K-mer size of 31 
# NOTE: Update the --tmpdir path to a valid temporary directory on your system 
# with sufficient storage space, or remove the flag to use default.
salmon index \
    -t $GENTROME \
    -i $OUT \
    --decoys $DECOYS \
    -k 31 \
    -p $SLURM_CPUS_PER_TASK \
    --tmpdir ${TEMP_WORK}/temp/salmon
