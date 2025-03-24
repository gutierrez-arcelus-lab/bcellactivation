library(tidyverse)
library(tximport)

if (!file.exists("data")) dir.create("./data")
if (!file.exists("plots")) dir.create("./plots")

# Sample meta data
samples_keep <- 
    read_tsv("../data/sample_decode.tsv")

meta_data <- 
    "../data/metadata_pooledreps.tsv" |>
    read_tsv(col_names = c("sample_id", "fq1", "fq2")) |>
    filter(sample_id %in% samples_keep$sample_name)

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

salmon_files <- 
    sprintf("../results/salmon_pooledreps/%s/quant.sf", meta_data$sample_id) |>
    setNames(meta_data$sample_id)

txi <- tximport(salmon_files, 
		type = "salmon", 
		tx2gene = select(tx_to_gene, tx_id, gene_id))

# Save data
write_rds(txi, "./data/gene_expression.rds")
