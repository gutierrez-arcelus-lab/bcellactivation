# README

## Step 0

Run setup.R to create a VCF with variants selected for ASE analysis.

## Step 1

Run subset\_mgb.slurm and merge\_mgb.slurm to obtain genotypes for 
the selected variants, for females of European ancestry in MGBB.

After running the pipeline, I realized that the region chr17:1-500,000
is missing from all chr17 VCFs. Regenerating the tabix index resolves
the issue. Therefore, I repeated the VCF subset process for chr17
copying the VCFs to $TEMP, generating a new index, and running the 
commands again.

This is how I created the array spec:

awk $2 == "chr17" ../../0-genotypes/array\_spec.txt > array\_spec\_chr17.txt

Then I run the subset and merge scripts again (saved as "\*_chr17.sh").

## Step 2

Run prep\_snp\_df.R to create a BED file with allele frequencies.

## Step 3

Run pileup.slurm to create QuASAR input files for every RNA-seq sample.

## Step 4

Run run\_quasar.R to run ASE analysis.

