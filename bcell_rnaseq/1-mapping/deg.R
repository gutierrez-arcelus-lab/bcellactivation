library(tidyverse)
library(DESeq2)
library(tximport)

# Salmon
salmon_files <- 
    list.files("./results", pattern = "quant.sf", recursive = TRUE, full.names = TRUE) 

names(salmon_files) <- 
    str_extract(salmon_files, "./results/(.+)/quant.sf", group = 1)

## select Day0 and DN2 samples
salmon_files <- salmon_files[grepl("unstday0|DN2", salmon_files)]

# Gene annotation
tx_annots <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.primary_assembly.annotation.gtf.gz" |>
    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc") |>
    filter(X3 == "transcript") |>
    mutate(transcript_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"), 
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	   gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_id = str_remove(gene_id, "\\.\\d+$")
	   ) |>
    select(transcript_id, gene_id, gene_name)

tx_to_gene <- select(tx_annots, TXNAME = transcript_id, GENEID = gene_id)

# Import expression data
txi <- tximport(salmon_files, type = "salmon", tx2gene = tx_to_gene)


# Run DEG analysis
sample_table <- 
    tibble(sample_id = colnames(txi$counts)) |>
    separate(sample_id, c("donor", "replic", "condition"), sep = "_", remove = FALSE) |>
    column_to_rownames("sample_id") |>
    mutate_all(factor)

dds <- DESeqDataSetFromTximport(txi, sample_table, ~ condition + donor)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "DN2", "unstday0"), alpha = 0.05)

as_tibble(res, rownames = "gene_id") |>
    filter(!is.na(padj)) |>
    left_join(distinct(tx_annots, gene_id, gene_name), join_by(gene_id)) |>
    select(gene_id, gene_name, everything()) |>
    arrange(padj) |>
    filter(padj <= 0.05)
