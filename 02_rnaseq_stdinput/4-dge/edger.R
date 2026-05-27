# ==============================================================================
# Description:  Performs Differential Gene Expression (DGE) analysis using edgeR. 
#               It imports transcript-level Salmon quantifications and aggregates 
#               them to the gene level using tximport. To allow for direct 
#		fold-change comparisons, this script uses the GENCODE v38 
#		annotation from the low-input pipeline.
# Input:        1. GENCODE v38 GTF (from low-input directory)
#               2. ./data/metadata.tsv (Sample sheet)
#               3. ./results/salmon/*/quant.sf (Salmon outputs)
# Output:       ./results/edger/results.tsv (Combined DGE summary table)
# ==============================================================================

library(tximport)
library(edgeR)
library(tidyverse)
library(glue)

# ------------------------------------------------------------------------------
# 1. Transcript-to-Gene Mapping (GENCODE v38)
# ------------------------------------------------------------------------------
# Import the exact same GENCODE v38 GTF used in the low-input pipeline
gtf <- 
    "../../01_rnaseq_lowinput/1_quantification/data/gencode.v38.primary_assembly.annotation.gtf.gz" |>
    rtracklayer::import(feature.type = "transcript")

# Extract the mapping table and strip the version numbers from the Ensembl IDs 
gene_tx <- 
    gtf |>
    as_tibble() |>
    select(transcript_id, gene_id, gene_name) |> 
    mutate(gene_id = str_remove(gene_id, "\\.\\d+"))

# ------------------------------------------------------------------------------
# 2. Metadata & File Parsing
# ------------------------------------------------------------------------------
sample_ids <- 
    "./data/metadata.tsv" |>
    read_tsv(col_types = "c--") |>
    pull(sample_id)

# Create the experimental design table.
# Splits the sample_id (e.g., BCR.24.donor1) into components and unites 
# stimulation and timepoint into a single 'group' variable.
sample_table <-
    tibble(sample_id = sample_ids) |>
    separate(sample_id, c("stim", "timep", "donor_id"), sep = "\\.", remove = FALSE) |>
    unite("group", c(stim, timep), sep = "_") |>
    column_to_rownames("sample_id")

# Create a named vector of paths to the Salmon quantification files
salmon_files <- 
    glue("./results/salmon/{sample_ids}/quant.sf") |>
    setNames(sample_ids)

# ------------------------------------------------------------------------------
# 3. Data Import via tximport
# ------------------------------------------------------------------------------
# Import Salmon transcript-level counts and aggregate to the gene level.
txi <- 
    tximport(salmon_files, 
	     type = "salmon", 
	     tx2gene = select(gene_tx, transcript_id, gene_id))

cts <- txi$counts
normMat <- txi$length

# ------------------------------------------------------------------------------
# 4. edgeR Offset Calculation
# ------------------------------------------------------------------------------
# To use Salmon counts in edgeR, we must calculate an offset matrix to account 
# for both transcript length changes across samples and library composition biases.

# Obtain per-observation scaling factors for length. 
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Compute effective library sizes from scaled counts using TMM normalization, 
# to account for composition biases between samples.
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combine effective library sizes with the length factors, 
# and calculate offsets for a log-link Generalized Linear Model (GLM).
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Create the DGEList object and apply the computed offset matrix
y <- DGEList(cts)
y <- scaleOffset(y, normMat)

# ------------------------------------------------------------------------------
# 5. edgeR Model Fitting & Filtering
# ------------------------------------------------------------------------------

# Define the design matrix (no intercept, allowing direct contrast definitions)
design <- model.matrix(~0 + group, data = sample_table)
colnames(design) <- sub("group", "", colnames(design))

# Filter out lowly expressed genes. 
# min.prop = 6/15 ensures a gene is kept if it has sufficient expression in a 
# minimum fraction of samples (e.g., equivalent to the smallest group size).
keep_y <- filterByExpr(y, design, large.n = 0, min.prop = 6/15)
y <- y[keep_y, ]

# Estimate dispersions and fit a Quasi-Likelihood (QL) negative binomial GLM.
# robust = TRUE protects empirical Bayes estimates against outlier genes.
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)

# ------------------------------------------------------------------------------
# 6. Differential Expression Testing
# ------------------------------------------------------------------------------

# Define the pairwise contrasts comparing each stimulation condition to baseline
my_contrasts <- 
    makeContrasts(TLR7 = TLR7_24 - unst_0,
		  BCR = BCR_24 - unst_0,
		  DN2 = DN2_24 - unst_0,
		  levels = design)

# Perform Quasi-Likelihood F-tests for each contrast and extract full tables
res_tlr7 <- 
    glmQLFTest(fit, contrast = my_contrasts[, "TLR7"]) |>
    topTags(n = Inf) |>
    {function(x) x$table}() |>
    rownames_to_column("gene_id") |>
    as_tibble()

res_bcr <- 
    glmQLFTest(fit, contrast = my_contrasts[, "BCR"]) |>
    topTags(n = Inf) |>
    {function(x) x$table}() |>
    rownames_to_column("gene_id") |>
    as_tibble()

res_dn2 <- 
    glmQLFTest(fit, contrast = my_contrasts[, "DN2"]) |>
    topTags(n = Inf) |>
    {function(x) x$table}() |>
    rownames_to_column("gene_id") |>
    as_tibble()

# ------------------------------------------------------------------------------
# 7. Output Compilation
# ------------------------------------------------------------------------------

# Bind all results together, merge in the gene names, and export the final table.
res <- 
    bind_rows("TLR7" = res_tlr7, "BCR" = res_bcr, "DN2" = res_dn2, .id = "stim") |>
    left_join(distinct(gene_tx, gene_id, gene_name)) |>
    select(stim, gene_id, gene_name, everything())

write_tsv(res, "./results/edger/results.tsv")

