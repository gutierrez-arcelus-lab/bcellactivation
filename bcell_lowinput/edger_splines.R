library(tximport)
library(edgeR)
library(splines)
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
meta_data <- 
    "./data/metadata_pooledreps.tsv" |>
    read_tsv(col_names = c("sample_id", "f1", "f2")) |>
    separate(sample_id, c("donor", "treat", "time"), sep = "_", remove = FALSE) |>
    filter(! treat %in% c("Unstim", "IL4")) |>
    mutate(treat = factor(treat, levels = c("CD40L", "TLR7", "TLR9", "BCR", "BCR-TLR7", "DN2")),
	   time = factor(time, levels = c("4hrs", "24hrs", "48hrs", "72hrs"))) |>
    arrange(treat, time)

salmon_files <- 
    sprintf("./results/salmon_pooledreps/%s/quant.sf", meta_data$sample_id) |>
    setNames(meta_data$sample_id)

# Import expression data
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

# Remove samples that failed
y <- y[, y$samples$lib.size > 2e6]

# Filtering using the design information
stims <- c("CD40L", "TLR9", "TLR7", "BCR", "BCR_TLR7", "DN2")

sample_table <- 
    y$samples |>
    as_tibble(rownames = "sample_id") |>
    separate(sample_id, c("donor", "treat", "time"), sep = "_", remove = FALSE) |>
    mutate(time = parse_number(time),
	   treat = recode(treat, "BCR-TLR7" = "BCR_TLR7"),
	   treat = factor(treat, levels = stims)) |>
    select(sample_id, donor, treat, time) |>
    arrange(treat, donor, time) |>
    column_to_rownames("sample_id")

y <- y[, rownames(sample_table)]

# Design
treatment <- sample_table$treat
X <- ns(parse_number(as.character(sample_table$time)), df = 3)

design <- model.matrix(~0 + treatment + treatment:X)
colnames(design) <- str_remove(colnames(design), "treatment")

keep_y <- filterByExpr(y, design)
y <- y[keep_y, ]


# Run edgeR
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)

# Test for all
const_mat <- matrix(0, nrow = 4, ncol = 6)

combinations <- combn(ncol(const_mat), 2)

gs <- colnames(design)

run_qlf <- function(i) {

    const_mat_i <- const_mat
    comb_i <- combinations[, i]
    const_mat_i[2:4, comb_i[1]] <- -1
    const_mat_i[2:4, comb_i[2]] <- 1
    contrasts_i <- const_mat_i |> t() |> as.integer()

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

write_tsv(qlf_all, "./results/edger/diff_expr_splines.tsv")