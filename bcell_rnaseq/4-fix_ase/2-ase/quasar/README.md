# README

## Step 0

Run setup.R to create a VCF with variants selected for ASE analysis.

## Step 1

Run subset_mgb.slurm and merge_mgb.slurm to obtain genotypes for 
the selected variants, for females of European ancestry in MGBB.

## Step 2

Run prep_snp_df.R to create a BED file with allele frequencies.

## Step 3

Run pileup.slurm to create QuASAR input files for every RNA-seq sample.

## Step 4

Run run_quasar.R to run ASE analysis.

