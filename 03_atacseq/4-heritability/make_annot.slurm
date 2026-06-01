#!/usr/bin/bash 

# ==============================================================================
# Description:  Creates binary SNP annotations and computes LD scores for our 
#               custom ATAC-seq peak sets. The Slurm array (1-110) iterates 
#               across 22 chromosomes for 5 distinct peak sets (4 stimulations + 
#               1 baseline control) in parallel.
# Input:        1. ./data/array_spec.tsv (Task grid from setup.R)
#               2. ./data/peak_coord.txt (All GRCh37 peak coordinates)
#               3. ./data/{SET}.txt (Condition-specific peak IDs)
#               4. 1000G Phase 3 BIM files and HapMap3 SNP list (institutional paths)
# Output:       Chromosome-specific .annot, .M, .M_5_50, and .l2.ldscore.gz 
#               files for each peak set in the ./data/ldscores/ directory.
# ==============================================================================

# ==============================================================================
# SLURM Configuration
# ==============================================================================

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=8gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-110
#SBATCH -J ldsc
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/ldsc-annot.%A-%a.out
#SBATCH -e logs/ldsc-annot.%A-%a.err

# ==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids (specific to our institutional cluster) and initialize Conda.
source /programs/biogrids.shrc
source ${HOME}/miniconda3/etc/profile.d/conda.sh

conda activate ldsc

cd $SLURM_SUBMIT_DIR

# ==============================================================================
# Input Parsing
# ==============================================================================
SPEC="./data/array_spec.tsv"

# Extract the target chromosome and peak set for this specific array task
CHROM=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $1 }' $SPEC)
SET=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'FNR==ARRID { print $2 }' $SPEC)

# Define paths to the required LDSC 1000G reference files
# NOTE FOR EXTERNAL USERS: Download these data from the LDSC repository
# See setup.R
BIM="./data/ldsc_data/v4/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHROM}.bim"
HAPMAP="./data/ldsc_data/v4/hm3_no_MHC.list.txt"

# Define custom peak data inputs and outputs
COORD="./data/peak_coord.txt"
PEAKSET="./data/${SET}.txt"
OUT="./data/ldscores/${SET}.${CHROM}"
ANNOT="${OUT}.annot"

# ==============================================================================
# 1. Create SNP Annotation
# ==============================================================================
# Use make_annot.py to identify which 1000G SNPs fall within our ATAC-seq peaks.
# --windowsize 500: Expands the peak by 500bp in both directions. This accounts 
#                   for slightly shifted causal variants that might lie just 
#                   outside the absolute center of the open chromatin region.
# NOTE: Requires a python environment with pybedtools installed.
python.pybedtools /lab-share/IM-Gutierrez-e2/Public/tools/ldsc/make_annot.py \
    --gene-set-file $PEAKSET \
    --gene-coord-file $COORD \
    --bimfile $BIM \
    --windowsize 500 \
    --annot-file $ANNOT

# ==============================================================================
# 2. Compute LD Scores
# ==============================================================================
# Use ldsc.py to calculate the actual LD scores for the newly annotated SNPs.
# --l2: Instructs the script to calculate LD scores.
# --print-snps: Restricts output strictly to HapMap3 SNPs (standard LDSC practice).
python2.7 /lab-share/IM-Gutierrez-e2/Public/tools/ldsc/ldsc.py \
    --l2 \
    --bfile ${BIM%.*} \
    --annot $ANNOT \
    --print-snps $HAPMAP \
    --ld-wind-cm 1 \
    --thin-annot \
    --out $OUT

# ==============================================================================
# 3. Rename Annotation Header
# ==============================================================================
# LDSC-SEG compares specific sets against the baseline. If multiple score files 
# have the exact same column header ("L2"), the downstream merge step will crash.
# This decompresses the score file, uses `sed` to rename the column header to 
# include the specific set name (e.g., "L2BCR"), and recompresses it.
zcat "${OUT}.l2.ldscore.gz" | sed "1,1s/L2/L2${SET}/" | gzip > "${OUT}.l2.tmp.gz"
mv "${OUT}.l2.tmp.gz" "${OUT}.l2.ldscore.gz" 

conda deactivate

