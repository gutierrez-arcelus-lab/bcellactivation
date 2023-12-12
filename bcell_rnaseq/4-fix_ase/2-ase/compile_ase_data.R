library(tidyverse)

if (!file.exists("plots")) dir.create("plots")

meta <- 
    "./array_spec.tsv" |>
    read_tsv(col_names = c("donor_id", "sample_id", "stim", "mgbid"), 
	     col_types = c(.default = "c")) |>
    unite("id", c("sample_id", "stim"), sep = "_", remove = FALSE) |>
    select(id, donor_id, sample_id, stim)
    
ase_df <- 
    sprintf("./results/%s.asereadcounter.txt", meta$id) |>
    setNames(meta$id) |>
    map_df(read_tsv, .id = "id") |>
    left_join(meta, by = "id") |>
    select(donor_id, sample_id, stim, everything()) |>
    select(-id) |>
    mutate(stim = recode(stim, "unstday0" = "Day 0"),
	   stim = factor(stim, levels = c("Day 0", "BCR", "TLR7", "DN2")))

sample_order <- ase_df |> 
    distinct(sample_id, stim) |> 
    count(sample_id, sort = TRUE) |>
    pull(sample_id)

ase_clean_df <- ase_df |>
    filter(totalCount >= 10,
	   otherBases/refCount < .1 & otherBases/altCount < .1,
	   (otherBases/(refCount + altCount)) < 0.05) |>
    select(sample_id, stim, var_id = variantID, refCount, altCount) |>
    mutate(sample_id = factor(sample_id, levels = sample_order))

annotations <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v39.primary_assembly.annotation.gtf") |>
    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc")

exon_annot <- annotations |>
    filter(X3 == "exon", X1 %in% paste0("chr", c(1:22, "X"))) |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(chr = X1, start = X4, end = X5, gene_id, gene_name)

gene_annot <- annotations |>
    filter(X3 == "gene", X1 %in% paste0("chr", c(1:22, "X"))) |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(chr = X1, start = X4, end = X5, gene_id, gene_name)

exonic_variants <- 
    ase_clean_df |>
    distinct(var_id) |>
    separate(var_id, c("chr", "pos", "ref", "alt"), 
	     sep = ":", convert = TRUE, remove = FALSE) |>
    inner_join(exon_annot, join_by(chr, between(pos, start, end))) |>
    distinct(var_id, gene_id, gene_name)

intronic_variants <- 
    ase_clean_df |>
    distinct(var_id) |>
    anti_join(exonic_variants, join_by(var_id)) |>
    separate(var_id, c("chr", "pos", "ref", "alt"), 
	     sep = ":", convert = TRUE, remove = FALSE) |>
    inner_join(gene_annot, join_by(chr, between(pos, start, end))) |>
    distinct(var_id, gene_id, gene_name)

variant_annot <- 
    bind_rows("exon" = exonic_variants, 
	      "intron" = intronic_variants,
	      .id = "annot") |>
    group_by(annot, var_id) |>
    summarise(gene_id = paste(gene_id, collapse = "/"),
	      gene_name = paste(gene_name, collapse = "/")) |>
    ungroup()

ase_annotated <- ase_clean_df |>
    left_join(variant_annot, join_by(var_id), relationship = "many-to-many") |>
    mutate(annot = ifelse(is.na(annot), "intergenic", annot)) |>
    select(sample_id, stim, var_id, annot, gene_id, gene_name, refCount, altCount)

write_tsv(ase_annotated, "./ase_data.tsv")
