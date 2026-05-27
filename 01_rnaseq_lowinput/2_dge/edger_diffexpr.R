# ==============================================================================
# Description:  Performs differential gene expression (DGE) analysis using edgeR.
#               Imports Salmon quantifications via tximport, calculates offsets 
#               to correct for sample-specific transcript length changes, and 
#               executes quasi-likelihood F-tests (QLF) across ALL pairwise 
#               combinations of experimental groups using mirai parallelization.
# Input:        1. gencode.v38.primary_assembly.annotation.gtf.gz (Annotations)
#               2. samples_pass.tsv (List of samples passing QC)
#               3. metadata.tsv (Sample sheet)
#               4. quant.sf files (Salmon quantifications)
# Output:       diff_expr_all_times_all_genes.tsv (Complete DGE results)
# ==============================================================================

# Running this with 36GB of RAM to handle the large dataframes
library(tximport)
library(edgeR)
library(tidyverse)
library(mirai)

# ------------------------------------------------------------------------------
# 1. Transcript to Gene Mapping
# ------------------------------------------------------------------------------
# Import the GTF using rtracklayer and extract only transcript features.
# This creates a tx2gene mapping table required by tximport to summarize 
# transcript-level counts to the gene level.
tx_to_gene <- 
    "../1_quantification/data/gencode.v38.primary_assembly.annotation.gtf.gz" |>
    rtracklayer::import(feature.type = "transcript") |>
    as_tibble() |>
    select(tx_id = transcript_id, gene_id, gene_name)

# ------------------------------------------------------------------------------
# 2. Experimental Design & Metadata Setup
# ------------------------------------------------------------------------------
stims <- c("Unstim", "IL4", "CD40L", "TLR9", "TLR7", "BCR", "BCR_TLR7", "DN2")

# Load the list of valid samples that passed the 2M read threshold in the QC step
samples_keep <- read_tsv("../0_qc/data/samples_pass.tsv")

# Parse the metadata to create a structured sample table
sample_table <- 
    "../data/metadata.tsv" |>
    read_tsv(col_names = c("sample_id", "f1", "f2")) |>
    select(sample_id) |>
    separate(sample_id, c("donor", "treat", "time"), sep = "_", remove = FALSE) |>
    mutate(treat = sub("-", "_", treat),
	   treat = factor(treat, levels = stims),
	   time = parse_number(time)) |>
    unite("group", c(treat, time), sep = ".", remove = FALSE) |>
    select(sample_id, group, donor, treat, time) |>
    arrange(treat, time, donor) |>
    mutate(group = fct_inorder(group)) |>
    filter(sample_id %in% samples_keep$sample_name) |>
    column_to_rownames("sample_id")

# ------------------------------------------------------------------------------
# 3. Import Data via tximport
# ------------------------------------------------------------------------------

# Define Salmon quantification files
salmon_files <- 
    sprintf("../1_quantification/results/%s/quant.sf", rownames(sample_table)) |>
    setNames(rownames(sample_table))

# Import transcript-level estimates and summarize to gene-level.
txi <- 
    tximport(salmon_files, 
	     type = "salmon", 
	     tx2gene = select(tx_to_gene, tx_id, gene_id))

cts <- txi$counts
normMat <- txi$length

# ------------------------------------------------------------------------------
# 4. edgeR Normalization with Length Offsets
# ------------------------------------------------------------------------------
# Standard edgeR assumes constant transcript length across samples. 
# Because alternative splicing can alter average gene lengths between conditions, 
# we use the tximport lengths to calculate an observation-specific offset matrix.

# 1. Obtain per-observation scaling factors for length 
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# 2. Compute effective library sizes from scaled counts, 
# to account for composition biases between samples.
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# 3. Combine effective library sizes with the length factors, 
# and calculate offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Create the edgeR DGEList object and inject the calculated offset matrix
y <- DGEList(cts)
y <- scaleOffset(y, normMat)

# Ensure sample table matches the columns of the DGEList
sample_table_filt <- 
    sample_table |>
    rownames_to_column("sample_id") |>
    filter(sample_id %in% colnames(y$counts)) |>
    column_to_rownames("sample_id")

y <- y[, rownames(sample_table_filt)]

# ------------------------------------------------------------------------------
# 5. Model Design & Filtering
# ------------------------------------------------------------------------------
# Specify the design matrix without an intercept (~0) to allow 
# comparisons between any two group means.
design <- model.matrix(~0 + group, data = sample_table_filt)
colnames(design) <- sub("group", "", colnames(design))

# Filter out lowly expressed genes based on the experimental design
keep_y <- filterByExpr(y, design, group = sample_table_filt$group)
y <- y[keep_y, ]

# Estimate dispersions and fit the Quasi-Likelihood (QL) negative binomial model
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)

# ------------------------------------------------------------------------------
# 6. Pairwise Testing Across All Combinations
# ------------------------------------------------------------------------------
gs <- colnames(design)

# Generate all possible pairwise contrasts between groups
contrast_strings <- 
    combn(gs, 2, FUN = function(x) paste(x[2], x[1], sep = "-"))

contrast_matrix <- 
    makeContrasts(contrasts = contrast_strings, levels = design)

# Use mirai to parallelize the QLF tests.
# This efficiently processes all pairwise comparisons while retaining all genes 
daemons(8)

dge_results_list <- 
    map(setNames(seq_along(contrast_strings), contrast_strings),
	in_parallel(function(i)
		    edgeR::glmQLFTest(fit, contrast = contrast_matrix[, i]) |>
		    edgeR::topTags(n = Inf) |>
		    {function(x) x$table}() |>
		    tibble::rownames_to_column("gene_id") |>
		    tibble::as_tibble(),
		    contrast_matrix = contrast_matrix, fit = fit),
    .progress = TRUE)

daemons(0)

# ------------------------------------------------------------------------------
# 7. Compile and Export Final Results
# ------------------------------------------------------------------------------
dge_results_df <- 
    dge_results_list |>
    bind_rows(.id = "contrast") |>
    left_join(distinct(tx_to_gene, gene_id, gene_name), join_by(gene_id)) |>
    select(contrast, gene_id, gene_name, everything())

write_tsv(dge_results_df, "./results/diff_expr_all_times_all_genes.tsv")
