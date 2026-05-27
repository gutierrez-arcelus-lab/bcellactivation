#!/usr/bin/bash

# ==============================================================================
# Description:  Generates a transcriptome FASTA file from the primary assembly. 
# Rationale:    At the time of the GENCODE v41 release, a pre-compiled 
#               transcriptome FASTA was only available for reference chromosomes. 
#               To ensure accurate mapping against the complete primary assembly, 
#               this script manually extracts the primary transcripts using RSEM.
# Input:        GRCh38.primary_assembly.genome.fa (Genome sequence)
#               gencode.v41.primary_assembly.annotation.gtf (Annotations)
# Output:       gencode.v41.primary_assembly.transcripts.fa
# Dependencies: Biogrids (RSEM v1.2.25)
# ==============================================================================

# ------------------------------------------------------------------------------
# SLURM DIRECTIVES
# ------------------------------------------------------------------------------

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --job-name=RSEM
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/RSEM.%j.out
#SBATCH -e logs/RSEM.%j.err

# ------------------------------------------------------------------------------
# ENVIRONMENT SETUP
# ------------------------------------------------------------------------------
# Note: Biogrids is our institutional software environment manager. 
# If running on a different HPC, replace these lines with your local equivalent.
source /programs/biogrids.shrc
export RSEM_X=1.2.25

# Move to the directory where the job was submitted from
cd "$SLURM_SUBMIT_DIR"

# ------------------------------------------------------------------------------
# INPUT / OUTPUT DEFINITIONS
# ------------------------------------------------------------------------------
# Define the reference primary assembly genome sequence
GENOME="./data/GRCh38.primary_assembly.genome.fa"

# Define the corresponding GENCODE v41 primary assembly annotations
GTF="./data/gencode.v41.primary_assembly.annotation.gtf"

# Define the output prefix for the RSEM reference builder
OUT="./data/gencode.v41.primary_assembly"

# ------------------------------------------------------------------------------
# EXECUTION & CLEANUP
# ------------------------------------------------------------------------------
# Run rsem-prepare-reference to parse the GTF and Genome to generate our 
# target transcript FASTA file (${OUT}.transcripts.fa).
rsem-prepare-reference --gtf "$GTF" "$GENOME" "$OUT"

# RSEM automatically generates several index and grouping files required for 
# its own mapping tools. Because we are only using RSEM to extract the 
# transcripts.fa file for downstream tools (e.g., Salmon), we remove these 
# unnecessary index files.
rm "$OUT.grp" "$OUT.ti" "$OUT.chrlist" "$OUT.seq" "$OUT.idx.fa" "$OUT.n2g.idx.fa"
