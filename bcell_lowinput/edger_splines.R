library(tximport)
library(edgeR)
library(splines)
library(tidyverse)

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

# Use Unstim 0hrs to create a base level for each stim
salmon_files <- 
    sprintf("./results/salmon_pooledreps/%s/quant.sf", meta_data$sample_id) |>
    setNames(meta_data$sample_id)

salmon_files_tmp <- 
    salmon_files |>
    enframe("sample_id", "file") |>
    separate(sample_id, c("donor", "treat", "time"), sep = "_", remove = FALSE) |>
    filter(treat != "IL4", !(treat == "Unstim" & time != "0hrs"))

salmon_files_baselevel <-
    salmon_files_tmp |>
    filter(treat == "Unstim", time == "0hrs") |>
    select(-treat) |>
    expand_grid(distinct(salmon_files_tmp, treat) |> filter(treat != "Unstim")) |>
    select(sample_id, donor, treat, time, file) |>
    unite("sample_id", c(donor, treat, time), sep = "_", remove = FALSE)

salmon_files_recoded <- 
    bind_rows(filter(salmon_files_tmp, treat != "Unstim"), salmon_files_baselevel) |>
    unite("group", c(treat, time), sep = "_", remove = FALSE) |>
    mutate(time = factor(time, levels = c("0hrs", "4hrs", "24hrs", "48hrs", "72hrs"))) |>
    arrange(treat, donor, time) |>
    select(sample_id, file) |>
    deframe()

# Import expression data
txi <- tximport(salmon_files_recoded, 
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
    mutate(donor = fct_inorder(donor),
	   time = fct_inorder(time),
	   treat = factor(treat)) |>
    select(sample_id, group, donor, treat, time) |>
    column_to_rownames("sample_id")

y <- y[, rownames(sample_table)]

# Splines
donor <- sample_table$donor
treatment <- sample_table$treat
X <- ns(parse_number(as.character(sample_table$time)), df = 3)

# Design
design <- model.matrix(~donor + treatment + treatment:X)
keep_y <- filterByExpr(y, design, group = sample_table$group)
y <- y[keep_y, ]

colnames(design)

# Run edgeR
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)




