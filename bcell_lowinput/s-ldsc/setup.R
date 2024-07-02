library(tidyverse)

gene_sets <- 
    "../wgcna/data/DN2_kme.tsv" |>
    read_tsv() |>
    select(-grey) |>
    pivot_longer(-gene_id, names_to = "module", values_to = "kme") |>
    group_by(gene_id) |>
    slice_max(kme) |>
    ungroup() |>
    arrange(module, desc(kme))

# Gene annotations
annot <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v38.primary_assembly.annotation.gtf.gz") |>
    read_tsv(comment = "#", col_names = FALSE) |>
    filter(X3 == "gene") |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(chr = X1, start = X4, end = X5, gene_id, gene_name, strand = X7)

# Create BED file
bed38 <- 
    gene_sets |>
    left_join(annot, join_by(gene_id)) |>
    filter(chr %in% paste0("chr", c(1:22))) |>
    group_by(module) |>
    top_n(250, kme) |>
    ungroup() |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    select(chr, start, end, gene_id, module) |>
    mutate(chr = factor(chr, levels = paste0("chr", 1:22)),
	   start = start - 1L) |>
    arrange(chr, start, end)

# lift over to hg19
bed38_file <- "./data/gene_sets/hg38.bed"
bed19_file <- "./data/gene_sets/hg19.bed"
fail_file <- "./data/gene_sets/failTolift.txt"
chain_file <- "/reference_databases/ReferenceGenome/liftover_chain/hg38/hg38ToHg19.over.chain.gz"

write_tsv(bed38, bed38_file, col_names = FALSE)

glue::glue("liftOver {bed38_file} {chain_file} {bed19_file} {fail_file}") |>
    system()

bed19 <- 
    read_tsv(bed19_file, col_names = c("chr", "start", "end", "gene_id", "module"))

out <- 
    bed19 |>
    select(-gene_id) |>
    group_by(chr, module) |>
    nest() |>
    ungroup() |>
    mutate(chr = str_remove(chr, "chr"),
	   f = glue::glue("./data/gene_sets/bed/{module}.{chr}.bed"))

# slurm array specification
out |>
    select(chr, module) |>
    mutate(chr = factor(chr, levels = 1:22),
	   module = factor(module)) |>
    complete(chr, module) |>
    write_tsv("./data/gene_sets/array_spec.tsv", col_names = FALSE)

# Save a BED file for each gene set
walk2(out$data, out$f, ~write_tsv(.x, .y, col_names = FALSE))

# Save Gene IDS for each module

gene_list <- 
    bed19 |>
    select(module, gene_id) |>
    {function(x) split(x, x$module)}() |>
    map(~pull(., gene_id))

walk2(gene_list, names(gene_list), 
     ~write_lines(.x, glue::glue("./data/gene_sets/genes/{.y}.txt")))

# Gene TSS file
strand_df <- 
    annot |>
    select(gene_id, strand) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

bed19 |>
    left_join(strand_df, join_by(gene_id)) |>
    mutate(TSS = case_when(strand == "+" ~ start,
			     strand == "-" ~ end),
	   START = TSS - 1L) |>
    select(GENE = gene_id, CHR = chr, START, END = TSS) |>
    write_tsv("./data/gene_sets/ENSG_coord.txt")
