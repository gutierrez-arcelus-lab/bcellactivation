#!/usr/bin/bash

#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=12gb
#SBATCH -t 08:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-4
#SBATCH -J homer
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/homer.%A-%a.out
#SBATCH -e logs/homer.%A-%a.err

# ==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids (Specific to our institutional cluster).
# External users: Ensure HOMER is installed and available in your $PATH.
source /programs/biogrids.shrc
export HOMER_X=4.11

cd "${SLURM_SUBMIT_DIR}"

# ==============================================================================
# Input Parsing
# ==============================================================================
# Extract the specific stimulation condition for this array task (e.g., "BCR")
STIM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' ./data/stims.txt )

# Define paths to the inputs required for HOMER
# Note: Ensure the reference genome FASTA is present or symlinked in this directory.
GENOME="./data/GRCh38.primary_assembly.genome.fa"
BED="./data/${STIM}.bed"
BG="./data/consensus.bed"

# Define output directory
OUT="./results/${STIM}"

# ==============================================================================
# Motif Enrichment Analysis
# ==============================================================================
# Run findMotifsGenome.pl:
# Positional arguments: 1) Target BED, 2) Genome FASTA, 3) Output directory
# -bg: The custom background BED file (consensus peaks) to prevent GC-content bias.
# -size 200: Extracts 200bp centered on the peak (optimal for narrow ATAC peaks).
# -mask: Ignores repeat-masked (lowercase) sequences in the reference genome.
# -mset vertebrates: Restricts the motif database strictly to vertebrate TFs.
# -p: Number of CPU threads to use for parallel processing.
findMotifsGenome.pl "$BED" "$GENOME" "$OUT" \
    -bg "$BG" \
    -size 200 \
    -mask \
    -mset vertebrates \
    -p "$SLURM_CPUS_PER_TASK"
