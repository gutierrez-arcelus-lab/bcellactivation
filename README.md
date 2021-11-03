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

<img src="./plots/pca_bcell_expression.png" width="2181" />

#### Fold change in comparison with “resting”

<img src="./plots/fc.png" width="2181" />

#### Subset of the plot above:

<img src="./plots/fc_subset.png" width="2181" />

### TO DO:

-   Run kallisto + sleuth once we have more samples, and call
    significant genes

## MGB Biobank analysis
