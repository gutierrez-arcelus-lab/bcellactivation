#!/usr/bin/bash 

# ==============================================================================
# Description:  Executes Partitioned Heritability analysis using LD Score 
#               Regression. It evaluates 37 distinct GWAS traits against the 
#               ATAC-seq peaks to determine if specific B cell states are 
#               enriched for disease heritability. It runs two models:
#               1. S-LDSC: Joint enrichment of all stimulated states.
#               2. LDSC-SEG: Specific enrichment of states vs. a baseline control.
# Input:        1. ./data/traits.txt (List of 37 GWAS traits)
#               2. Pre-calculated LD scores (./data/ldscores/*.ldscore.gz)
#               3. ./data/peaksets.ldcts (Continuous trait score configuration)
#               4. Baseline 1000G reference weights and frequencies.
# Output:       S-LDSC (`sldsc_*`) and LDSC-SEG (`ldscseg_*`) results files in 
#               the ./results/ directory.
# ==============================================================================

# ==============================================================================
# SLURM Configuration
# ==============================================================================

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=16gb
#SBATCH -t 24:00:00
#SBATCH -p bch-compute
#SBATCH -a 1-37
#SBATCH -J ldsc
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/ldsc.%A-%a.out
#SBATCH -e logs/ldsc.%A-%a.err

# ==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids (specific to our institutional cluster) and initialize Conda.
# LDSC strictly requires a Python 2.7 environment to run.
source /programs/biogrids.shrc
source ${HOME}/miniconda3/etc/profile.d/conda.sh

conda activate ldsc

cd $SLURM_SUBMIT_DIR

# ==============================================================================
# Input Parsing (GWAS Traits)
# ==============================================================================
METADATA="./data/traits.txt"

# Extract the LDSC sub-directory and the trait identifier based on the array ID.
STATSDIR=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" -F '\t' 'FNR==ARRID { print $1 }' $METADATA ) 
TRAIT=$( awk -v ARRID="$SLURM_ARRAY_TASK_ID" -F '\t' 'FNR==ARRID { print $3 }' $METADATA )

# ==============================================================================
# Define Reference Paths
# ==============================================================================
# NOTE FOR EXTERNAL USERS: 
# You will need to download the "1000G Phase 3" baseline models, weights, and 
# frequencies from the broadinstitute/ldsc GitHub/Zenodo repository and update 
# these paths accordingly.
# see setup.R
LDSCDIR="./data/ldsc_data/v4" 
GWAS="${LDSCDIR}/${STATSDIR}/${TRAIT}.sumstats.gz" 
WEIGHTS="${LDSCDIR}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
FREQ="${LDSCDIR}/1000G_Phase3_frq/1000G.EUR.QC." 
BASELINE="${LDSCDIR}/1000G_Phase3_baseline_v1.2_ldscores/baseline."

# Define prefixes for our custom calculated ATAC-seq LD scores
IL4="./data/ldscores/IL4."
TLR="./data/ldscores/TLR7."
BCR="./data/ldscores/BCR."
DN2="./data/ldscores/DN2."
CTS="./data/peaksets.ldcts"

# ==============================================================================
# 1. Stratified LD Score Regression (S-LDSC)
# ==============================================================================
# Runs a joint model analyzing all four stimulation states simultaneously against 
# the baseline model.
# --overlap-annot: Required because our custom peak sets might overlap with 
#                  annotations in the baseline model.
python2.7 /lab-share/IM-Gutierrez-e2/Public/tools/ldsc/ldsc.py \
    --h2 $GWAS \
    --ref-ld-chr $BASELINE,$IL4,$TLR,$BCR,$DN2 \
    --w-ld-chr $WEIGHTS \
    --frqfile-chr $FREQ \
    --overlap-annot --print-coefficients \
    --out ./results/sldsc_${TRAIT}

# ==============================================================================
# 2. Cell-Type-Specific LD Score Regression (LDSC-SEG)
# ==============================================================================
# Runs the specific expression/accessibility model to test if our highly 
# accessible peaks (vs. the custom control background defined in the CTS file) 
# drive disease heritability.
python2.7 /lab-share/IM-Gutierrez-e2/Public/tools/ldsc/ldsc.py \
    --h2-cts $GWAS \
    --ref-ld-chr $BASELINE \
    --ref-ld-chr-cts $CTS \
    --w-ld-chr $WEIGHTS \
    --out ./results/ldscseg_${TRAIT}

conda deactivate
