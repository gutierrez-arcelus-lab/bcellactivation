# Workflow

## 0. Setup

First, we need a specificiation table to parallelize jobs. This table is create by the script `make_array.sh`.

Also, we pre-select female individuals. First we run `compute_chrx_het.slurm` to compute individual-level heterozygosities. After visual inspection with the `../plot.R` script, we write the IDs to `./results/females_BATCH.txt` files.


## 1. Process VCFs

### 1.1. Filter VCFs from MGB and 1000 Genomes

Filtering is carried out by the scripts `filter_mgb_vcf.slurm` and `filter_1000G_vcf.slurm`.

### 1.2. Merge VCFs

Then we merge the VCFs from 1000G and MGB for each chromosome, for the set of variants present in all VCFs (intersect), is run LD pruning with plink. This step is carried out by the script `merge_mgb_1000G.slurm`.

### 1.3. Concatenate

The script `concat_vcf.slurm` concatenates VCFs for each chromosome into a single final VCF, which is saved in the `results` directory.

### 1.4. PCA

Submit the script `run_pca.slurm` to run PCA with plink.

### 1.5. Select MGB individuals with high European ancestry

Run the R script `select_eur_individuals.R`


## 2. Process SLE risk variants

We take the supplementary table from Langefeld et al. (2017) with 3 tiers with significant variants for different ancestries (directory `sle_variants`).


- With the R script `parse_sle_variants.R`, we select all variants at FDR < 5% for Europeans;
- Lift the original coordinates in hg19 to hg38 with the script `liftsle.sh`;
- Extract the genotypes for the selected variants from each chromosome and merge into a single VCF file with `extract_sle_vcf.sh`.


## Bonus: ADMIXTURE

- Run unsupervised analysis on the 1000 Genomes data with `admixture_1000G_vc.slurm`;
- Project MGB samples onto 1000G reference with `admixture_projection.sh`.
