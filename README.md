README
================

## B-Cells: Bulk RNA-seq

### Input data:

-   1 individual, 4 timepoints:
    -   Resting for 16 hours
    -   IgG stim for 24 hours
    -   IgG stim for 72 hours
    -   RSQ stim for 24 hours
    -   RSQ stim for 72 hours

Fastq files located on directory:
/lab-share/IM-Gutierrez-e2/Public/B\_cells/bulkTCpilot\_1/34.198.31.178/210618\_MG8989\_fastq

### Methods

-   Expression levels were estimated with Salmon (script
    `./salmon_quant.slurm`)

### Results

#### Overview of expression levels

<img src="./plots/transcript_biotypes.png" width="2181" />

#### PCA of expression levels

<img src="./plots/pca_bcell_expression.png" width="2275" />

#### Fold change in comparison with “resting”

<img src="./plots/fc.png" width="2181" />

#### Subset of the plot above:

<img src="./plots/fc_subset.png" width="2181" />

#### Correlations with resting state:

<img src="./plots/scatter_resting_conditions.png" width="2228" />

### TO DO:

-   Run kallisto + sleuth once we have more samples, and call
    significant genes

## MGB Biobank analysis

### Input data

-   4921 individuals;
-   \~79M variants.

### Methods

The procedure below is carried out by the `./process_vcf.slurm` script:

-   VCF processing:
    -   remove variants with any missing genotypes;
    -   select only biallelic SNPs with MAF &gt;= 0.1;
    -   remove variants with pvalue &lt; 0.05/79M for HWE test.
-   1000 Genomes data:
    -   low coverage and exome data realigned to GRCh38 (not NYGC
        version);
    -   same processing as above, except that I did not remove variants
        failing HWE.
-   Match MGB and 1000G VCFs:
    -   remove A/T and C/G genotypes due to potential strand ambiguity;
    -   remove duplicate entries (these can be multiallelic variants);
    -   select variants with the same position and alleles in both
        datasets;
    -   filter both datasets for the common set of variants;
    -   merge VCFs and run LD pruning for r2 &gt; 0.2.
-   PCA:
    -   QTLtools pca

### Results

### TO DO:

-   select SLE variants from Langefeld et al. (2017)
-   select individuals from MGB biobank
