#!/usr/bin/bash

# ==============================================================================
# Description:  Aligns adapter-trimmed paired-end RNA-seq reads to the human 
#               genome using STAR in 2-pass mode to improve novel junction 
#               discovery. Extracts high-confidence splice junctions from 
#               uniquely mapped reads using Regtools for downstream splicing analysis.
# Input:        `metadata.tsv` (Format: sample_id, fastq_R1, fastq_R2). 
#               Note: Multiple fastq files for a single sample should 
#               be provided as a comma-separated string in the fastq columns.
# Output:       1. Coordinate-sorted BAM files
#               2. Extracted splice junction files (.junc)
# ==============================================================================

# ==============================================================================
# SLURM Configuration
# ==============================================================================

#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=48gb
#SBATCH --time=48:00:00 
#SBATCH -p bch-compute 
#SBATCH --array=1-74
#SBATCH --job-name=STAR
#SBATCH --mail-user=vitor.aguiar@childrens.harvard.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/STAR-%A-%a.out
#SBATCH -e logs/STAR-%A-%a.err

# ==============================================================================
# Environment Setup
# ==============================================================================
# Load BioGrids (Specific to our institutional cluster).
# External users: Ensure STAR, Regtools, and SAMtools are in your $PATH.
source /programs/biogrids.shrc
export STAR_X=2.7.11a
export REGTOOLS_X=0.5.1
export SAMTOOLS_X=1.13

cd $SLURM_SUBMIT_DIR

# ==============================================================================
# Input Parsing
# ==============================================================================
META="./data/metadata.tsv"

# Parse the sample ID and comma-separated FASTQ lists.
# `awk 'FNR>1'` skips the header row.
SAMPLEID=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $1 }' )
FQ1=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $2 }' ) 
FQ2=$( awk 'FNR>1' $META | awk -v ARRID="$SLURM_ARRAY_TASK_ID" 'NR==ARRID { print $3 }' ) 

# Define paths
INDEX="./index"
OUT="./bam/${SAMPLEID}" 

# ==============================================================================
# 1. STAR Mapping
# ==============================================================================
# Run STAR alignment:
# --readFilesIn: Natively accepts comma-separated lists to merge lanes on the fly.
# --twopassMode Basic: Runs a 1st pass to discover novel junctions, then maps 
#                      reads in a 2nd pass using those novel junctions.
# --outFilterMismatchNoverReadLmax 0.04: Restricts mismatches to max 4% of read length.
STAR --runMode alignReads \
    --runThreadN "$SLURM_CPUS_PER_TASK" \
    --genomeDir "$INDEX" \
    --readFilesIn "$FQ1" "$FQ2" \
    --readFilesCommand zcat  \
    --outFilterMismatchNoverReadLmax 0.04 \
    --twopassMode Basic \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${OUT}_"
   
# ==============================================================================
# 2. Splice Junction Extraction (for downstream Splicing Analysis)
# ==============================================================================
BAM="${OUT}_Aligned.sortedByCoord.out.bam"

# Index the raw output BAM
samtools index -@ "$SLURM_CPUS_PER_TASK" "$BAM"

# Filter for strictly uniquely mapping reads on primary autosomes + chrX.
# STAR uses MAPQ=255 to indicate uniquely mapped reads.
samtools view -@ "$SLURM_CPUS_PER_TASK" -b -h -q 255 "$BAM" `(seq 1 22; echo X) | sed 's/^/chr/'` > "${OUT}.uniq.bam"
samtools index -@ "$SLURM_CPUS_PER_TASK" "${OUT}.uniq.bam"

# Extract junctions using Regtools:
# -s 0: Unstranded data
# -a 8: Minimum anchor length of 8bp on each side of the junction
# -m 50 -M 500000: Minimum/Maximum intron size
regtools junctions extract -s 0 -a 8 -m 50 -M 500000 "${OUT}.uniq.bam" -o "${OUT}.junc"

# ==============================================================================
# 3. Cleanup
# ==============================================================================
# Remove the intermediate unique BAM and temporary STAR directories
rm "${OUT}.uniq.bam"*
rm -r "${OUT}__"*
