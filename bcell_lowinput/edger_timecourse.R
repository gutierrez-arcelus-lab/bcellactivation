library(tximport)
library(edgeR)
library(splines)
library(tidyverse)

# Functions
run_edger <- function(stim, dge) {

    # subset for stim of interest
    # using Unstim_0hrs as starting point for all stims
    pattern <- sprintf("Unstim_0hrs|_%s_", stim)
    dge_s <- dge[, grepl(pattern, colnames(dge))] 

    # get timepoint of each sample
    time_vector <- 
	str_split(colnames(dge_s), "_") |>
	map_chr(3) |>
	str_remove("hrs$") |>
	as.integer()
   
    dge_s$samples$group <- factor(time_vector, levels = sort(unique(time_vector)))

    # Create design matrix
    X <- ns(time_vector, df = 3)
    design <- model.matrix(~X)

    # Estimate dispersion
    dge_s <- estimateDisp(dge_s, design)
   
    # Time course trend analysis
    fit <- dge_s |>
	glmQLFit(design, robust = TRUE) |>
	glmQLFTest(coef = 2:4)

    # Obtain results
    res_df <- topTags(fit, n = Inf) |>
	as.data.frame() |>
	as_tibble(rownames = "gene_id")

    # Obtain observed and fitted CPM values
    obs_cpm <- 
	cpm(dge_s, offset = dge_s$offset, log = FALSE) |>
	as_tibble(rownames = "gene_id") |>
	pivot_longer(-gene_id, names_to = "sample_id", values_to = "obs_cpm")
   
    # this fit cpm is not on the same scale as the observed cpm
    # try to fix it
    fit_cpm <- 
	cpm(fit, log = FALSE, lib.size = exp(getOffset(dge_s))) |>
	`colnames<-`(colnames(dge_s$counts)) |>
	as_tibble(rownames = "gene_id") |>
	pivot_longer(-gene_id, names_to = "sample_id", values_to = "fit_cpm")

    cpm_df <- left_join(obs_cpm, fit_cpm, join_by(gene_id, sample_id))

    # Return results
    list("results" = res_df, "cpm" = cpm_df)
}

# Sample meta data
meta_data <- read_tsv("./data/sample_decode.tsv") |>
    separate(sample_name, c("name", "stim", "time"), sep = "_", remove = FALSE) |>
    mutate(time = factor(time, levels = str_sort(unique(time), numeric = TRUE))) |>
    arrange(stim, time, name)

# Transcript to Gene map
tx_to_gene <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v38.primary_assembly.annotation.gtf") |>
    read_tsv(comment = "#", col_names = FALSE) |>
    filter(X3 == "transcript") |>
    mutate(tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	   gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(tx_id, gene_id, gene_name)

# Import expression data
salmon_files <- 
    sprintf("./results/salmon/%s/quant.sf", meta_data$sample_id) |>
    setNames(meta_data$sample_id)

txi <- tximport(salmon_files, 
		type = "salmon", 
		tx2gene = select(tx_to_gene, tx_id, gene_id))

cts <- txi$counts
normMat <- txi$length

# Obtaining per-observation scaling factors for length, 
# adjusted to avoid changing the magnitude of the counts.
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

# Filtering using the design information
sample_table <- meta_data |>
    unite("group", c(stim, time), sep = ".") |>
    select(sample_id, group) |>
    column_to_rownames("sample_id")

design <- model.matrix(~group, data = sample_table)
keep_y <- filterByExpr(y, design)
y <- y[keep_y, ]

# For samples with technical replicates, 
# keep only the one with highest library size.
# I could pooled technical replicates with sumTechReps,
# but need to verify its behavior when combining offsets.
#pooled <- sumTechReps(y, ID = meta_data$sample_name)
reps_to_rm <- y$samples |>
    as_tibble(rownames = "sample_id") |>
    left_join(meta_data) |>
    add_count(sample_name) |>
    filter(n > 1) |>
    group_by(sample_name) |>
    slice(-which.max(lib.size)) |>
    ungroup() |>
    pull(sample_id)

y$counts <- y$counts[, !colnames(y$counts) %in% reps_to_rm]
y$samples <- y$samples[colnames(y$counts), ]
y$offset <- y$offset[, colnames(y$counts)]

# Rename samples
sample_names <- 
    tibble(sample_id = colnames(y$counts)) |>
    left_join(meta_data) |>
    pull(sample_name)

colnames(y$counts) <- rownames(y$samples) <- colnames(y$offset) <- sample_names

# Run edgeR
stims <- meta_data |>
    distinct(stim) |>
    filter(! stim %in% c("Unstim", "IL4")) |>
    pull() |>
    {function(x) setNames(x, x)}()

results_list <- map(stims, run_edger, dge = y)

# Save results
results_df <- map_df(results_list, "results", .id = "stim") |>
    left_join(distinct(tx_to_gene, gene_id, gene_name), join_by(gene_id)) |>
    select(stim, gene_id, gene_name, everything())

cpm_df <- map_df(results_list, "cpm", .id = "stim") |>
    left_join(distinct(tx_to_gene, gene_id, gene_name), join_by(gene_id)) |>
    select(stim, gene_id, gene_name, everything())

write_tsv(cpm_df, "./results/edger/cpm.tsv")
write_tsv(results_df, "./results/edger/results.tsv")
