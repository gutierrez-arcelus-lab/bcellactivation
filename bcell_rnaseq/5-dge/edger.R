library(tximport)
library(edgeR)
library(tidyverse)
library(glue)

# Transcript annotations
gtf <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v41.primary_assembly.annotation.gtf.gz" |>
    rtracklayer::import(feature.type = "transcript")

gene_tx <- 
    gtf |>
    as.data.frame() |>
    as_tibble() |>
    select(transcript_id, gene_id, gene_name) |> 
    mutate(gene_id = str_remove(gene_id, "\\.\\d+"))

sample_ids <- 
    "../4-splicing/data/metadata_qced_pooled.tsv" |>
    read_tsv(col_types = "c--") |>
    pull(sample_id)

sample_table <-
    tibble(sample_id = sample_ids) |>
    separate(sample_id, c("stim", "timep", "donor_id"), sep = "\\.", remove = FALSE) |>
    unite("group", c(stim, timep), sep = "_") |>
    column_to_rownames("sample_id")

salmon_files <- 
    glue("../4-splicing/salmon_quant_pooled/{sample_ids}/quant.sf") |>
    setNames(sample_ids)

txi <- 
    tximport(salmon_files, 
	     type = "salmon", 
	     tx2gene = select(gene_tx, transcript_id, gene_id))

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

design <- model.matrix(~0 + group, data = sample_table)

colnames(design) <- sub("group", "", colnames(design))

keep_y <- filterByExpr(y, design, group = sample_table$group)
y <- y[keep_y, ]

# Run edgeR
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design, robust = TRUE)

my_contrasts <- 
    makeContrasts(TLR7 = TLR7_24 - unst_0,
		  BCR = BCR_24 - unst_0,
		  DN2 = DN2_24 - unst_0,
		  levels = design)

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

res <- bind_rows("TLR7" = res_tlr7, "BCR" = res_bcr, "DN2" = res_dn2, .id = "stim")

write_tsv(res, "./results.tsv")

