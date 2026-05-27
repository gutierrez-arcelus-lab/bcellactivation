#!/usr/bin/bash

# ==============================================================================
# Description:  Generates a STAR genome index using the GRCh38 primary assembly
#               and GENCODE v41 annotations. This index is required for downstream
#               splicing-aware read mapping.
# Input:        1. GRCh38.primary_assembly.genome.fa (Primary assembly FASTA)
#               2. gencode.v41.primary_assembly.annotation.gtf (Annotation)
# Output:       Compiled STAR genome index in the `./index` directory.
# ==============================================================================

# ==============================================================================
# SLURM Configuration
# ==============================================================================

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --time=24:00:00
#SBATCH --mem=48gb
#SBATCH -p bch-compute
#SBATCH --job-name=STARidx
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBACTH --mail-type=END,FAIL
#SBATCH -o logs/STAR-index.o%A-%a.out
#SBATCH -e logs/STAR-index.o%A-%a.err

# ==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids (Specific to our institutional cluster).
# External users: Ensure STAR is installed and available in your $PATH.
source /programs/biogrids.shrc
export STAR_X=2.7.11a

cd $SLURM_SUBMIT_DIR

# ==============================================================================
# Variable Definition
# ==============================================================================
GENOME="./data/GRCh38.primary_assembly.genome.fa"
GTF="./data/gencode.v41.primary_assembly.annotation.gtf"
OH=149
OUT="./index"

mkdir -p "$OUT"

# ==============================================================================
# Generate STAR Index
# ==============================================================================
STAR --runMode genomeGenerate \
    --runThreadN "$SLURM_CPUS_PER_TASK" \
    --genomeDir "$OUT" \
    --genomeFastaFiles "$GENOME" \
    --sjdbGTFfile "$GTF" \
    --sjdbOverhang "$OH"
