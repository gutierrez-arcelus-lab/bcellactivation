library(edgeR)
library(splines)
library(tidyverse)

# Functions
run_edger <- function(count_matrix) {

    # get timepoint of each sample
    time_vector <- sub("^[^_]+_(\\d+)hrs$", "\\1", colnames(count_matrix)) |>
	as.integer()
    
    # Create the DGEList object, filter and normalized data
    dge <- DGEList(counts = count_matrix, group = time_vector)
    keep_y <- filterByExpr(dge)
    dge <- dge[keep_y, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)

    # Create design matrix
    X <- ns(time_vector, df = 3)
    design <- model.matrix(~X)

    # Estimate dispersion
    dge <- estimateDisp(dge, design)
   
    # Time course trend analysis
    fit <- dge |>
	glmQLFit(design, robust = TRUE) |>
	glmQLFTest(coef = 2:4)

    # Obtain results
    res_df <- topTags(fit, n = Inf) |>
	as.data.frame() |>
	as_tibble(rownames = "gene_id")

    # Obtain observed and fitted CPM values
    obs_cpm <- 
	cpm(dge, log = FALSE, prior.count = fit$prior.count) |>
	as_tibble(rownames = "gene_id") |>
	pivot_longer(-gene_id, names_to = "sample_id", values_to = "obs_cpm")
    
    fit_cpm <- 
	cpm(fit, log = FALSE) |>
	`colnames<-`(colnames(dge$counts)) |>
	as_tibble(rownames = "gene_id") |>
	pivot_longer(-gene_id, names_to = "sample_id", values_to = "fit_cpm")

    cpm_df <- left_join(obs_cpm, fit_cpm, join_by(gene_id, sample_id))

    # Return results
    list("results" = res_df, "cpm" = cpm_df)

}



# Expression data
expr_df <- read_rds("./data/expression.rds") |>
    mutate(gene_id = sub("\\.\\d+$", "", gene_id))

# Gene info
genes_df <- distinct(expr_df, gene_id, gene_name)

# Sample meta data
meta_data <- read_tsv("./data/sample_decode.tsv") |>
    separate(sample_name, c("name", "stim", "time"), sep = "_", remove = FALSE) |>
    mutate(time = factor(time, levels = str_sort(unique(time), numeric = TRUE))) |>
    arrange(stim, time, name)

# Count matrix
counts <- expr_df |>
    select(-tpm, -gene_name) |>
    pivot_wider(names_from = id, values_from = count) |>
    select(gene_id, all_of(meta_data$sample_id)) |>
    column_to_rownames("gene_id") |>
    data.matrix()

# Pool technical replicates
pooled <- sumTechReps(counts, ID = meta_data$sample_name)

# Use Unstim 0hr as baseline for all stims
matrix_0hrs <- 
    pooled[ , grep("Unstim_0hrs$", colnames(pooled), value = TRUE)] |>
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "sample_name") |>
    left_join(distinct(meta_data, sample_name, name, time), join_by(sample_name)) |>
    select(-sample_name) |>
    pivot_wider(names_from = time, values_from = value)

matrix_stim <- 
    pooled[ , !grepl("Unstim_0hrs$", colnames(pooled))] |>
    as_tibble(rownames = "gene_id") |>
    pivot_longer(-gene_id, names_to = "sample_name") |>
    left_join(distinct(meta_data, sample_name, name, stim, time), join_by(sample_name)) |>
    select(-sample_name) |>
    pivot_wider(names_from = time, values_from = value)
    
# Make final expression matrices
count_matrix_list <- 
    left_join(matrix_stim, matrix_0hrs, join_by(gene_id, name)) |>
    filter(! stim %in% c("Unstim", "IL4")) |>
    select(gene_id, name, stim, `0hrs`, everything()) |>
    pivot_longer(-c(gene_id, name, stim), names_to = "time", values_to = "count") |>
    drop_na(count) |>
    unite("sample_id", c(name, time), sep = "_") |>
    {function(x) split(x, x$stim)}() |>
    map(~select(., -stim) |>
	pivot_wider(names_from = sample_id, values_from = count) |>
	column_to_rownames("gene_id") |>
	data.matrix())

# Run edgeR
results_list <- map(count_matrix_list, run_edger)

# Save results
results_df <- map_df(results_list, "results", .id = "stim")
cpm_df <- map_df(results_list, "cpm", .id = "stim")

write_tsv(cpm_df, "./results/edger/cpm.tsv")
write_tsv(results_df, "./results/edger/results.tsv")
