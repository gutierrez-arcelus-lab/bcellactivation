library(tidyverse)

gene_sets <- 
    "../wgcna/data/DN2_kme.tsv" |>
    read_tsv() |>
    select(-grey) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    pivot_longer(-gene_id, names_to = "module", values_to = "kme") |>
    #group_by(gene_id) |>
    #slice_max(kme) |>
    #ungroup() |>
    arrange(module, desc(kme))

# Gene annotations
annot19 <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v19.chr_patch_hapl_scaff.annotation.gtf") |>
    read_tsv(comment = "#", col_names = FALSE) |>
    filter(X1 %in% paste0("chr", 1:22), X3 == "gene") |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	   gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    select(chr = X1, start = X4, end = X5, gene_id, gene_name, strand = X7)

annot38 <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v38.primary_assembly.annotation.gtf.gz") |>
    read_tsv(comment = "#", col_names = FALSE) |>
    filter(X1 %in% paste0("chr", 1:22), X3 == "gene") |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	   gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    select(chr = X1, start = X4, end = X5, gene_id, gene_name, strand = X7)

# Run liftOver for genes not annotated in hg19
bed38 <- 
    gene_sets |>
    distinct(gene_id) |>
    left_join(annot19, join_by(gene_id)) |>
    filter(is.na(start)) |>
    select(gene_id) |>
    inner_join(annot38, join_by(gene_id)) |>
    select(chr, start, end, gene_id) |>
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

bedlift <- 
    read_tsv(bed19_file, col_names = c("chr", "start", "end", "gene_id"))

# Merge and create final BED
bed_out <- 
    gene_sets |>
    distinct(gene_id) |>
    inner_join(annot19, join_by(gene_id)) |>
    select(chr, start, end, gene_id) |>
    bind_rows(bedlift) |>
    mutate(chr = factor(chr, levels = paste0("chr", 1:22))) |>
    arrange(chr, start, end)

# Gene TSS file
strand_df <- 
    select(annot38, gene_id, strand)

bed_out |>
    left_join(strand_df, join_by(gene_id)) |>
    mutate(TSS = case_when(strand == "+" ~ start,
			     strand == "-" ~ end),
	   START = TSS - 1L) |>
    select(GENE = gene_id, CHR = chr, START, END = TSS) |>
    write_tsv("./data/gene_sets/ENSG_coord.txt")

# slurm array specification
module_order <- c("turquoise", "blue", "brown", "yellow", "green", "red", "black")

expand_grid(chr = 1:22, 
	    module = c(unique(gene_sets$module), "control")) |>
    mutate(module = factor(module, levels = c(module_order, "control"))) |>
    arrange(chr, module) |>
    write_tsv("./data/gene_sets/array_spec.tsv", col_names = FALSE)

# Write CTS file
tibble(module = module_order) |>
    mutate(annot = glue::glue("data/gene_sets/ldscores/{module}."),
	   control = "data/gene_sets/ldscores/control.") |>
    unite("ldscores", c(annot, control), sep = ",") |>
    write_tsv("./data/gene_sets/module.ldcts", col_names = FALSE)

# Write gene IDS for each module
write_lines(bed_out$gene_id, "./data/gene_sets/genes/control.txt")

gene_list <- 
    gene_sets |>
    group_by(module) |>
    top_n(500, kme) |>
    ungroup() |>
    left_join(bed_out, join_by(gene_id)) |>
    arrange(chr, start, end) |>
    select(module, gene_id) |>
    {function(x) split(x, x$module)}() |>
    map(~pull(., gene_id))

walk2(gene_list, names(gene_list), 
     ~write_lines(.x, glue::glue("./data/gene_sets/genes/{.y}.txt")))
