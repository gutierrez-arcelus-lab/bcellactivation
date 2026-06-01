#!/usr/bin/bash

# ==============================================================================
# Description:  Demultiplexes pooled CITE-seq/scRNA-seq samples back to their 
#               original human donors using natural genetic variation (SNPs). 
#               It uses demuxlet to match single-cell BAM reads against a 
#               donor genotype VCF file.
# Input:        1. ../../data/cellranger/{library_ID}/.../bam (Cell Ranger BAM)
#               2. ../../data/cellranger/{library_ID}/.../barcodes.tsv.gz
#               3. ./data/citeseq_donors.vcf.gz (Donor genotypes)
# Output:       Demultiplexing assignments and doublet probabilities in 
#               the ./results/{LIB}/ directory.
# ==============================================================================

# ==============================================================================
# SLURM Configuration
# ==============================================================================

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH -p bch-compute
#SBATCH --array=1-3
#SBATCH --job-name=demuxlet
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/demuxlet_%A_%a.out
#SBATCH -e logs/demuxlet_%A_%a.err

# ==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids (Specific to our institutional cluster).
# External users: Ensure demuxlet is installed and available in your $PATH.
source /programs/biogrids.shrc

cd $SLURM_SUBMIT_DIR

# ==============================================================================
# Array Task Routing
# ==============================================================================
# Map the Slurm array ID to the specific library name.
case "${SLURM_ARRAY_TASK_ID}" in
    1) LIB="1984" ;;
    2) LIB="1988" ;;
    3) LIB="1990" ;; 
    *) echo "Invalid array task ID" && exit 1 ;;
esac

# ==============================================================================
# Define Relative Paths
# ==============================================================================

# External users: Ensure dbGaP downloads are extracted matching this structure,
# or adjust paths.
BAM="../../data/cellranger/${LIB}/BRI_${LIB}_possorted_genome_bam.bam"
BARCODE="../../data/cellranger/${LIB}/filtered_feature_bc_matrix/barcodes.tsv.gz"
VCF="./data/citeseq_donors.vcf.gz"

# Define and create output directory
RES="${SLURM_SUBMIT_DIR}/results/${LIB}"
mkdir -p "$RES"

# ==============================================================================
# Execute Demuxlet
# ==============================================================================
# Run demuxlet:
# --field GT: Uses the Genotype field from the VCF.
# --min-BQ 20: Ignores reads with a base quality score below 20 at the SNP position.
# --alpha: Defines the expected proportion of SNP errors (0) and doublet mixing (0.5).
demuxlet --sam "$BAM" \
    --vcf "$VCF" --field GT \
    --group-list "$BARCODE" \
    --min-BQ 20 \
    --alpha 0 --alpha 0.5 \
    --out "${RES}/results"
