library(tidyverse)
library(glue)

# Data directories
if (!file.exists("data")) dir.create("data")
if (!file.exists("data/ldscores")) dir.create("data/ldscores")
if (!file.exists("results")) dir.create("results")

# Import ASE data
gencode <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v45.primary_assembly.annotation.gtf.gz") |>
    vroom::vroom(comment = "#", col_names = FALSE, col_types = "ccciicccc") |>
    filter(X3 == "gene", X1 %in% paste0("chr", 1:22)) |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+")) |>
    select(chr = X1, start = X4, end = X5, gene_id)

ase_data <- 
    read_tsv("../3-ase/ase_data.tsv") |>
    filter(!is.na(gene_id))

### GLM for stim-dependent ASE
ase_glm <- read_tsv("../3-ase/results_glm/glm_res_df.tsv")


# Run liftOver to convert GRCh38 coordinates to GRCh37
bed38_file <- "./data/hg38.bed"
bed19_file <- "./data/hg19.bed"
fail_file <- "./data/failTolift.txt"
chain_file <- "/reference_databases/ReferenceGenome/liftover_chain/hg38/hg38ToHg19.over.chain.gz"

#bed_38 <- 
#    ase_data |>
#    distinct(gene_id) |>
#    rowid_to_column("idx") |>
#    separate_rows(gene_id, sep = "/") |>
#    inner_join(gencode, join_by(gene_id)) |>
#    group_by(idx) |>
#    slice(1) |>
#    ungroup() |>
#    mutate(chr = factor(chr, levels = paste0("chr", 1:22)),
#	   gene_id = str_remove(gene_id, "\\.\\d+$")) |>
#    arrange(chr, start, end) |>
#    select(chr, start, end, gene_id) |>
#    distinct()
#
bed_38 <- 
    ase_glm |>
    distinct(variant_id) |>
    left_join(distinct(ase_data, variant_id, gene_id), join_by(variant_id)) |>
    distinct(gene_id) |>
    rowid_to_column("idx") |>
    separate_rows(gene_id, sep = "/") |>
    inner_join(gencode, join_by(gene_id)) |>
    group_by(idx) |>
    slice(1) |>
    ungroup() |>
    mutate(chr = factor(chr, levels = paste0("chr", 1:22)),
	   gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    arrange(chr, start, end) |>
    select(chr, start, end, gene_id) |>
    distinct()

write_tsv(bed_38, bed38_file, col_names = FALSE)

glue("liftOver {bed38_file} {chain_file} {bed19_file} {fail_file}") |>
    system()

bedlift <- read_tsv(bed19_file, col_names = c("chr", "start", "end", "gene_id"))

# Make gene-coord-file for make_annot
bedlift |>
    select(GENE = gene_id, CHR = chr, START = start, END = end) |>
    write_tsv("./data/gene_coord.txt")

# Write IDs for control
write_lines(bedlift$gene_id, "./data/control.txt")

# Make gene sets
#gene_sets <- 
#    ase_data |>
#    filter(q_value <= 0.01) |>
#    #filter(p_bonferroni < 0.05) |>
#    distinct(stim, gene_id) |>
#    separate_rows(gene_id, sep = "/") |>
#    mutate(stim = str_remove(stim, " "),
#	   stim = factor(stim, levels = c("Day0", "TLR7", "BCR", "DN2")),
#	   gene_id = str_remove(gene_id, "\\.\\d+$")) |>
#    filter(gene_id %in% bedlift$gene_id) |>
#    {function(x) split(x, x$stim)}() |> 
#    map(~pull(., gene_id))
#
#walk2(gene_sets, names(gene_sets), 
#     ~write_lines(.x, glue("./data/{.y}.txt")))
#
#

#glm_top <- 
#    ase_glm |>
#    group_by(donor_id, replic, stim) |>
#    mutate(pct = ntile(abs(z), 100)) |>
#    ungroup() |> 
#    filter(pct >= 90) |> 
#    arrange(desc(abs(z)))

glm_top <- 
    ase_glm |>
    filter(grepl("^Day 0", stim)) |>
    mutate(stim = str_remove(stim, "Day 0-")) |>
    left_join(distinct(ase_data, variant_id, gene_id)) |>
    group_by(stim, gene_id) |>
    slice_max(abs(z)) |>
    ungroup() |>
    mutate(pct = ntile(abs(z), 100)) |>
    filter(pct >= 90) |>
    separate_rows(gene_id, sep = "/") |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    filter(gene_id %in% bedlift$gene_id) |>
    distinct(stim, gene_id)

gene_sets <- 
    glm_top |>
    mutate(stim = factor(stim, levels = c("TLR7", "BCR", "DN2"))) |>
    {function(x) split(x, x$stim)}() |> 
    map(~pull(., gene_id))

walk2(gene_sets, names(gene_sets), 
     ~write_lines(.x, glue("./data/{.y}.txt")))

# GWAS traits
read_tsv("../../atacseq/ldsc/data/traits.txt", col_names = FALSE) |>
    filter(X2 != "Cancer") |>
    write_tsv("./data/traits.tsv", col_names = FALSE)

# slurm array specification
stim_order <- names(gene_sets)

expand_grid(chr = 1:22, 
	    set = c(stim_order, "control")) |>
    mutate(set = fct_inorder(set)) |>
    arrange(chr, set) |>
    write_tsv("./data/array_spec.tsv", col_names = FALSE)

# Write CTS file
tibble(set = stim_order) |>
    mutate(annot = glue("data/ldscores/{set}."),
	   control = "data/ldscores/control.") |>
    unite("ldscores", c(annot, control), sep = ",") |>
    write_tsv("./data/genesets.ldcts", col_names = FALSE)

