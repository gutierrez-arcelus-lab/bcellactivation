library(tidyverse)
library(glue)
library(qvalue)

# Annotations
annotations <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v45.primary_assembly.annotation.gtf.gz") |>
    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc")

exon_annot <- 
    annotations |>
    filter(X3 == "exon", X1 %in% paste0("chr", c(1:22, "X"))) |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(chr = X1, start = X4, end = X5, gene_id, gene_name) |>
    distinct()

gene_annot <- 
    annotations |>
    filter(X3 == "gene", X1 %in% paste0("chr", c(1:22, "X"))) |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
    select(chr = X1, start = X4, end = X5, gene_id, gene_name)

# Static ASE
meta <- 
    "../2-mbv/matching_results.tsv" |>
    read_tsv() |>
    separate(rnaseq_sample_id, c("rnaseq_donor_id", "rep_id"), sep = "_") |>
    unite("sample_id", c(vcf_donor_id, rep_id), sep = "_") |>
    select(sample_id, stim = rnaseq_stim)
    
ase_df <-
    meta |>
    mutate(f = glue("./results/{sample_id}_{stim}_ase.txt"),
	   data = map(f, read_tsv)) |>
    mutate(stim = recode(stim, "unstday0" = "Day 0"),
	   stim = factor(stim, levels = c("Day 0", "BCR", "TLR7", "DN2"))) |>
    select(-f) |>
    unnest(cols = data)

ase_clean_df <- 
    ase_df |>
    filter(totalCount >= 20,
	   otherBases/refCount < .1,
	   otherBases/altCount < .1,
	   otherBases/totalCount < 0.05) |>
    select(sample_id, stim, variantID, refCount, altCount) |>
    janitor::clean_names()


exonic_variants <- 
    ase_clean_df |>
    distinct(variant_id) |>
    separate(variant_id, c("chr", "pos", "ref", "alt"), 
	     sep = ":", convert = TRUE, remove = FALSE) |>
    inner_join(exon_annot, join_by(chr, between(pos, start, end))) |>
    distinct(variant_id, gene_id, gene_name)

intronic_variants <- 
    ase_clean_df |>
    distinct(variant_id) |>
    anti_join(exonic_variants, join_by(variant_id)) |>
    separate(variant_id, c("chr", "pos", "ref", "alt"), 
	     sep = ":", convert = TRUE, remove = FALSE) |>
    inner_join(gene_annot, join_by(chr, between(pos, start, end))) |>
    distinct(variant_id, gene_id, gene_name)

variant_annot <- 
    bind_rows("exon" = exonic_variants, 
	      "intron" = intronic_variants,
	      .id = "annot") |>
    group_by(annot, variant_id) |>
    summarise(gene_id = paste(gene_id, collapse = "/"),
	      gene_name = paste(gene_name, collapse = "/")) |>
    ungroup()

ase_annotated <- 
    ase_clean_df |>
    left_join(variant_annot, join_by(variant_id), relationship = "many-to-many") |>
    mutate(annot = ifelse(is.na(annot), "intergenic", annot))

ase_res <- 
    ase_annotated |>
    mutate(p_value = map2_dbl(ref_count, ref_count + alt_count, 
			      ~binom.test(.x, .y, p = .5, alternative = "two.sided")$p.value),
	   q_value = qvalue(p_value)$qvalues)

write_tsv(ase_res, "./ase_data.tsv")


# Dynamic ASE
ase_data_dyn <- read_tsv("./results_glm/glm_res_df.tsv", col_types = "cfccdd") 

exonic_variants_dyn <- 
    ase_data_dyn |>
    distinct(variant_id) |>
    separate(variant_id, c("chr", "pos", "ref", "alt"), 
	     sep = ":", convert = TRUE, remove = FALSE) |>
    inner_join(exon_annot, join_by(chr, between(pos, start, end))) |>
    distinct(variant_id, gene_id, gene_name)

intronic_variants_dyn <- 
    ase_data_dyn |>
    distinct(variant_id) |>
    anti_join(exonic_variants_dyn, join_by(variant_id)) |>
    separate(variant_id, c("chr", "pos", "ref", "alt"), 
	     sep = ":", convert = TRUE, remove = FALSE) |>
    inner_join(gene_annot, join_by(chr, between(pos, start, end))) |>
    distinct(variant_id, gene_id, gene_name)

variant_annot_dyn <- 
    bind_rows("exon" = exonic_variants_dyn, 
	      "intron" = intronic_variants_dyn,
	      .id = "annot") |>
    group_by(annot, variant_id) |>
    summarise(gene_id = paste(gene_id, collapse = "/"),
	      gene_name = paste(gene_name, collapse = "/")) |>
    ungroup()

ase_dyn_annotated <- 
    ase_data_dyn |>
    left_join(variant_annot_dyn, join_by(variant_id), relationship = "many-to-many") |>
    mutate(annot = ifelse(is.na(annot), "intergenic", annot))

write_tsv(ase_dyn_annotated, "./results_glm/glm_res_df_annotated.tsv")

#ase_dyn_annotated <- read_tsv("./results_glm/glm_res_df_annotated.tsv")

# For LDSC
tested_genes <- 
    ase_dyn_annotated |>
    filter(grepl("^Day 0-", stim)) |>
    filter(!is.na(gene_id)) |>
    separate_rows(gene_id, gene_name, sep = "/") |>
    distinct(stim, gene_id, gene_name) |>
    arrange(stim, gene_name, gene_id)

signif_genes <-
    ase_dyn_annotated |>
    filter(grepl("^Day 0-", stim)) |>
    filter(!is.na(gene_id)) |>
    filter(p.adjust(p, method = "fdr") <= 0.1) |>
    separate_rows(gene_id, gene_name, sep = "/") |>
    distinct(stim, gene_id, gene_name) |>
    arrange(stim, gene_name, gene_id)

