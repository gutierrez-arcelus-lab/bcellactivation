# ==============================================================================
# Description:  Performs within-group time-course differential expression analysis.
#               Uses natural cubic splines to model non-linear temporal gene 
#               expression changes for each stimulation condition relative to a 
#               shared baseline (Unstim_0hrs).
# Input:        1. gencode.v38.primary_assembly.annotation.gtf.gz (Annotations)
#               2. metadata.tsv (Sample sheet)
#               3. quant.sf files (Salmon quantifications)
# Output:       1. results.tsv (Significant time-course trends per stimulation)
#               2. cpm.tsv (Observed CPM/logCPM values for plotting)
# ==============================================================================

library(tximport)
library(edgeR)
library(splines)
library(tidyverse)

# ------------------------------------------------------------------------------
# 1. Define Time-Course Analysis Function
# ------------------------------------------------------------------------------
# This function isolates a single stimulation, builds a spline-based model 
# across its time points, and tests for significant temporal dynamics.
run_edger <- function(stim, dge) {

    # Subset the DGE object for the current stimulation of interest, 
    # anchoring it to the shared Unstim_0hrs baseline.
    pattern <- sprintf("Unstim_0hrs|_%s_", stim)
    dge_s <- dge[, grepl(pattern, colnames(dge))] 

    # Extract the numeric timepoint from the sample names
    time_vector <- 
	str_split(colnames(dge_s), "_") |>
	map_chr(3) |>
	str_remove("hrs$") |>
	as.integer()
   
    # Ensure groups are factorized in chronological order
    dge_s$samples$group <- factor(time_vector, levels = sort(unique(time_vector)))

    # Create a dynamic design matrix using natural cubic splines (ns).
    # Degrees of freedom (df) are scaled based on available time points to 
    # capture non-linear trends.
    degrees <- 
	case_when(n_distinct(time_vector) <= 3 ~ 2L,
		  n_distinct(time_vector) > 3 ~ 3L,
		  TRUE ~ NA_integer_)

    X <- ns(time_vector, df = degrees)
    design <- model.matrix(~X)

    # Estimate sample dispersions
    dge_s <- estimateDisp(dge_s, design)
   
    # Time course trend analysis:
    # By testing all spline coefficients simultaneously (coef = 2:ncol(design)), 
    # we test the null hypothesis that the expression curve is completely flat 
    # over time. A significant p-value indicates the gene changed over time.
    fit <- 
	dge_s |>
	glmQLFit(design, robust = TRUE) |>
	glmQLFTest(coef = 2:ncol(design))

    # Extract results
    res_df <- 
	topTags(fit, n = Inf) |>
	as.data.frame() |>
	as_tibble(rownames = "gene_id")

    # Extract observed logCPM and CPM values to facilitate downstream visualization
    obs_logcpm <- 
	cpm(dge_s, offset = dge_s$offset, log = TRUE) |>
	as_tibble(rownames = "gene_id") |>
	pivot_longer(-gene_id, names_to = "sample_id", values_to = "obs_logcpm")
    
    obs_cpm <- 
	cpm(dge_s, offset = dge_s$offset, log = FALSE) |>
	as_tibble(rownames = "gene_id") |>
	pivot_longer(-gene_id, names_to = "sample_id", values_to = "obs_cpm")
   
    # Return both the statistics and expression estimates
    list("results" = res_df, "cpm" = left_join(obs_cpm, obs_logcpm, join_by(gene_id, sample_id)))
}

# ------------------------------------------------------------------------------
# 2. Transcript to Gene Mapping
# ------------------------------------------------------------------------------
tx_to_gene <- 
    "../1_quantification/data/gencode.v38.primary_assembly.annotation.gtf.gz" |>
    rtracklayer::import(feature.type = "transcript") |>
    as_tibble() |>
    select(tx_id = transcript_id, gene_id, gene_name)

# ------------------------------------------------------------------------------
# 3. Import Data via tximport
# ------------------------------------------------------------------------------
meta_data <- 
    "../data/metadata.tsv" |>
    read_tsv(col_names = c("sample_id", "f1", "f2"))

salmon_files <- 
    sprintf("../1_quantification/results/%s/quant.sf", meta_data$sample_id) |>
    setNames(meta_data$sample_id)

txi <- 
    tximport(salmon_files, 
	     type = "salmon", 
	     tx2gene = select(tx_to_gene, tx_id, gene_id))

cts <- txi$counts
normMat <- txi$length

# ------------------------------------------------------------------------------
# 4. edgeR Normalization with Length Offsets
# ------------------------------------------------------------------------------
# Obtaining per-observation scaling factors for length
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Computing effective library sizes from scaled counts, 
# to account for composition biases between samples.
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, 
# and calculating offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts)
y <- scaleOffset(y, normMat)

# ------------------------------------------------------------------------------
# 5. Global Filtering
# ------------------------------------------------------------------------------
# Remove samples that failed the initial read depth QC
y <- y[, y$samples$lib.size > 2e6]

# Filter out lowly expressed genes globally using a broad design matrix
sample_table <- 
    y$samples |>
    as_tibble(rownames = "sample_id") |>
    extract(sample_id, "group", "[^_]+_([^_]+_[^_]+)", remove = FALSE) |>
    select(sample_id, group) |>
    column_to_rownames("sample_id")

design <- model.matrix(~group, data = sample_table)
keep_y <- filterByExpr(y, design)
y <- y[keep_y, ]

# ------------------------------------------------------------------------------
# 6. Execute Time-Course Analysis per Stimulation
# ------------------------------------------------------------------------------
# Identify all unique stimulations (excluding the Unstim baseline)
stims <- 
    meta_data |>
    extract(sample_id, "stim", "[^_]+_([^_]+)_") |>
    distinct(stim) |>
    filter(! stim %in% c("Unstim")) |>
    pull() |>
    {function(x) setNames(x, x)}()

# Map the run_edger function across all stimulations to generate spline models
results_list <- map(stims, run_edger, dge = y)

# ------------------------------------------------------------------------------
# 7. Compile and Export Final Results
# ------------------------------------------------------------------------------
# Combine statistics
results_df <- 
    map_df(results_list, "results", .id = "stim") |>
    left_join(distinct(tx_to_gene, gene_id, gene_name), join_by(gene_id)) |>
    select(stim, gene_id, gene_name, everything())

# Combine observed CPMs
cpm_df <- 
    map_df(results_list, "cpm", .id = "stim") |>
    left_join(distinct(tx_to_gene, gene_id, gene_name), join_by(gene_id)) |>
    select(stim, gene_id, gene_name, everything())

write_tsv(cpm_df, "./results/edger/cpm.tsv")
write_tsv(results_df, "./results/edger/results.tsv")
