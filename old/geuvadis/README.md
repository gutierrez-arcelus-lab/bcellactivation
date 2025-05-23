README
================

## Trans-eQTL analysis

### 0. Define the individuals

Run the `write_sample_ids.R` script to extract individuals with EBV
data, and write a table with known covariates.

### 1. Filter the 1000 Genomes VCFs

Execute sbatch `filter_1000G_vcf.slurm`

in order to extract the genotype data for the individuals of interest,
and filter variants according the MAF and other parameters.

### 2. Make VCF with all chromosomes

Execute the `concat_vcf.sh` script to concatenate VCFs for each
chromosome into a single VCF.

### 3. PCA

Execute `sbatch pca.slurm` to perform PCA on genotype data.

### 4. Covariates

Execute the `make_covars.R` script to make a covariates matrix to be
used with QTLtools.

### 5. Prepare phenotype data

Execute the `prepare_qtl_pheno.R` to put the phenotype data in the
format which QTLtools expects.

### 6. Normalize expression

Execute the `qtltools_correct.sh` script to correct the expression data
for known covariates, and normalize to std normal.

### 7. Run trans eQTL nominal pass

Execute `sbatch qtltools_trans_nominal.slurm` to run the trans eQTL
analysis.

### 8. Run permutation

Execute `sbatch qtltools_trans_permutation.slurm` to run the
permutation.

### 9. Compute FDR

Execute the `qltools_trans_fdr.sh` script to compute the adjusted
p-values.

## Results

<img src="./plots/ebv_genes.png" width="3000" />

<img src="./plots/qqplot.png" width="1500" />

<img src="./plots/trans_manhattan.png" width="3059" />

### LMP-2B expression \~ chr18:71,147,350 SNP genotypes

<img src="./plots/trans_qtl.png" width="1800" />

### Correlation of expression human vs EBV

Black points = Significant after Bonferroni correction at alpha = 0.01

<img src="./plots/ebv_volcano.png" width="1200" />

<img src="./plots/ebv_heatmap.png" width="3600" />

<img src="./plots/ebv_pairs_corr.png" width="1800" />
