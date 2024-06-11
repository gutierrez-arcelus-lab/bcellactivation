library(tximport)
library(edgeR)
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
    degrees <- case_when(n_distinct(time_vector) <= 3 ~ 2L,
			 n_distinct(time_vector) > 3 ~ 3L,
			 TRUE ~ NA_integer_)

    X <- ns(time_vector, df = degrees)
    design <- model.matrix(~X)

    # Estimate dispersion
    dge_s <- estimateDisp(dge_s, design)
   
    # Time course trend analysis
    fit <- dge_s |>
	glmQLFit(design, robust = TRUE) |>
	glmQLFTest(coef = 2:ncol(design))

    # Obtain results
    res_df <- topTags(fit, n = Inf) |>
	as.data.frame() |>
	as_tibble(rownames = "gene_id")

    # Obtain observed and fitted CPM values
    obs_logcpm <- 
	cpm(dge_s, offset = dge_s$offset, log = TRUE) |>
	as_tibble(rownames = "gene_id") |>
	pivot_longer(-gene_id, names_to = "sample_id", values_to = "obs_logcpm")
    
    obs_cpm <- 
	cpm(dge_s, offset = dge_s$offset, log = FALSE) |>
	as_tibble(rownames = "gene_id") |>
	pivot_longer(-gene_id, names_to = "sample_id", values_to = "obs_cpm")
   
    # Return results
    list("results" = res_df, "cpm" = left_join(obs_cpm, obs_logcpm, join_by(gene_id, sample_id)))
}


# Transcript to Gene map
tx_to_gene <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v38.primary_assembly.annotation.gtf.gz") |>
    read_tsv(comment = "#", col_names = FALSE) |>
    filter(X3 == "transcript") |>
    mutate(tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
	   gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(tx_id, gene_id, gene_name)

# Import expression data
meta_data <- read_tsv("./data/metadata_pooledreps.tsv", col_names = c("sample_id", "f1", "f2"))

salmon_files <- 
    sprintf("./results/salmon_pooledreps/%s/quant.sf", meta_data$sample_id) |>
    setNames(meta_data$sample_id)

txi <- 
    tximport(salmon_files, 
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

# Remove samples that failed
y <- y[, y$samples$lib.size > 2e6]

# Filtering using the design information
sample_table <- 
    y$samples |>
    as_tibble(rownames = "sample_id") |>
    separate(sample_id, c("donor", "treat", "time"), sep = "_", remove = FALSE) |>
    unite("group", c(treat, time), sep = "_", remove = FALSE) |>
    filter(! treat %in% c("Unstim", "IL4")) |>
    mutate(time = factor(time, levels = c("4hrs", "24hrs", "48hrs", "72hrs")),
	   treat = factor(treat)) |>
    select(sample_id, group, donor, treat, time) |>
    arrange(treat, donor, time) |>
    column_to_rownames("sample_id")

y <- y[, rownames(sample_table)]

# Specify design
design <- model.matrix(~donor + treat + treat:time, data = sample_table)
keep_y <- filterByExpr(y, design, group = sample_table$group)
y <- y[keep_y, ]

# Run edgeR
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)

colnames(design)

qlf <- glmQLFTest(fit, coef = 27)
topTags(qlf)

