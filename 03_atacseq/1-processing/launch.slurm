#!/usr/bin/bash

# ==============================================================================
# Description:  Script to launch the nf-core/atacseq pipeline. 
#               This job acts as the Nextflow manager/orchestrator; it does not 
#               perform the heavy computation itself, but rather submits individual 
#               pipeline steps (mapping, peak calling, etc.) to the cluster.
# Input:        1. samplesheet.csv (nf-core formatted sample sheet)
#               2. e2.config (Custom cluster profile/configuration)
#               3. ./data/encode/blacklist/hg38-blacklist.v3.bed (Genome blacklist)
# Output:       Processed BAMs, narrow peaks, and MultiQC reports in `./results/`
# ==============================================================================

# ==============================================================================
# SLURM Configuration (Manager Node)
# ==============================================================================

#SBATCH --nodes=1                
#SBATCH --cpus-per-task=1        
#SBATCH --mem=4G                 
#SBATCH --time=72:00:00
#SBATCH --partition=bch-compute
#SBATCH --job-name=atac
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --output=logs/atac_nfcore.%j.out
#SBATCH --error=logs/atac_nfcore.%j.err

# ==============================================================================
# Environment Setup
# ==============================================================================
# Clear any conflicting modules and load Singularity.
# The nf-core pipeline relies on Singularity to pull and run software containers.
module purge
module load singularity

# Define cache and temporary directories for Nextflow and Singularity.
# NOTE FOR EXTERNAL USERS: You MUST update these to valid paths on your cluster
# with substantial storage capacity, as intermediate files can be very large.
export NXF_SINGULARITY_CACHEDIR="${TEMP_WORK}/temp/singularity/nxf_cache"
export SINGULARITY_TMPDIR="${TEMP_WORK}/temp/singularity/tmp"
export NXF_TEMP="${TEMP_WORK}/temp/nextflow/tmp"
export NXF_WORK="${TEMP_WORK}/temp/nextflow/work"

# Move to the directory where the job was submitted from
cd "$SLURM_SUBMIT_DIR"

# ==============================================================================
# Execute nf-core/atacseq
# ==============================================================================
# Launch Nextflow from a local Conda environment. 
# NOTE FOR EXTERNAL USERS: Adjust the `conda run` path to point to your own 
# Nextflow installation, or simply use `nextflow run ...` if it is in your $PATH.
# 
# Pipeline parameters:
# -r 2.0: Locks the pipeline to version 2.0.
# -c e2.config: Applies our institutional cluster execution settings.
# --read_length 150: Optimizes MACS2 peak calling for 150bp reads.
# --narrow_peak: ATAC-seq primarily yields narrow regions of open chromatin.
# --mito_name chrM: Essential for calculating and filtering mitochondrial reads.
# --blacklist: Filters out artifact-prone genomic regions (e.g., centromeres).
# -resume: Allows the pipeline to restart from the last successful step if interrupted.

conda run -p "${HOME}/miniconda3/envs/nf-core" nextflow run nf-core/atacseq -r 2.0 \
    -c "./data/e2.config" \
    --input "./data/samplesheet.csv" \
    --read_length 150 \
    --narrow_peak \
    --mito_name chrM \
    --blacklist "./data/encode/blacklist/hg38-blacklist.v3.bed" \
    --outdir "./results" \
    --multiqc_title B_cell_atac_nfcore \
    -resume

