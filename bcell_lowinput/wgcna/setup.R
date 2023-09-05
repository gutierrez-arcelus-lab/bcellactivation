library(tidyverse)
library(tximport)

# Sample meta data
meta_data <- 
    "../data/metadata_pooledreps.tsv" |>
    read_tsv(col_names = c("sample_id", "fq1", "fq2"))

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


# Expression data
salmon_files <- 
    sprintf("../results/salmon_pooledreps/%s/quant.sf", meta_data$sample_id) |>
    setNames(meta_data$sample_id)

txi <- tximport(salmon_files, 
		type = "salmon", 
		tx2gene = select(tx_to_gene, tx_id, gene_id))

# Save data
if (!file.exists("data")) dir.create("./data")

write_rds(txi, "./data/gene_expression.rds")

# Lupus-associated genes
sle_genes <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Khunsriraksakul/Khunsriraksakul_SuppData2.xlsx" |>
    readxl::read_excel(skip = 1) |>
    distinct(`Mapped Gene`) |>
    pull(1)

write_lines(sle_genes, "./data/lupus_genes.txt")
