library(tximport)
library(edgeR)
library(tidyverse)
library(furrr)

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
stims <- c("Unstim", "IL4", "CD40L", "TLR9", "TLR7", "BCR", "BCR_TLR7", "DN2")

samples_keep <- read_tsv("./data/sample_decode.tsv")

sample_table <- 
    "./data/metadata_pooledreps.tsv" |>
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

salmon_files <- 
    sprintf("./results/salmon_pooledreps/%s/quant.sf", rownames(sample_table)) |>
    setNames(rownames(sample_table))

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

sample_table_filt <- 
    sample_table |>
    rownames_to_column("sample_id") |>
    filter(sample_id %in% colnames(y$counts)) |>
    column_to_rownames("sample_id")

y <- y[, rownames(sample_table_filt)]

# Specify design
design <- model.matrix(~0 + group, data = sample_table_filt)

colnames(design) <- sub("group", "", colnames(design))

colnames(design)

keep_y <- filterByExpr(y, design, group = sample_table_filt$group)
y <- y[keep_y, ]

# Run edgeR
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)

# Test for all
gs <- colnames(design)

combinations <- combn(seq_len(length(gs)), 2)

run_qlf <- function(i) {
    
    comb_i <- combinations[ , i]
    contrasts_i <- rep(0, length(gs)) 
    contrasts_i[comb_i[1]] <- -1
    contrasts_i[comb_i[2]] <- 1

    qlf <- glmQLFTest(fit, contrast = contrasts_i)

    res <- 
	topTags(qlf, n = Inf, p.value = 0.05) |> 
	as.data.frame() |>
	as_tibble(rownames = "gene_id")

    if (nrow(res) > 0) {

	res |>
	    mutate(group1 = gs[comb_i[1]], group2 = gs[comb_i[2]]) |>
	    select(group1, group2, everything())

    } else {
	
	tibble(group1 = gs[comb_i[1]], group2 = gs[comb_i[2]])
    }

}

plan(multisession, workers = availableCores())

qlf_all <- 
    future_map_dfr(1:ncol(combinations), run_qlf) |>
    left_join(distinct(tx_to_gene, gene_id, gene_name), join_by(gene_id)) |>
    select(group1, group2, gene_id, gene_name, everything())

write_tsv(qlf_all, "./results/edger/diff_expr_all_times.tsv")
