library(tidyverse)
library(glue)
library(qvalue)
library(furrr)

plan(multisession, workers = availableCores())

# Annotations
annotations <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v45.primary_assembly.annotation.gtf.gz") |>
    vroom::vroom(comment = "#", col_names = FALSE, col_types = "ccciicccc")

exon_annot <- 
    annotations |>
    filter(X3 == "exon", X1 %in% paste0("chr", c(1:22, "X"))) |>
    group_split(X1) |>
    future_map_dfr(~mutate(., gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
			   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
		   select(chr = X1, start = X4, end = X5, gene_id, gene_name) |>
		   distinct())

gene_annot <- 
    annotations |>
    filter(X3 == "gene", X1 %in% paste0("chr", c(1:22, "X"))) |>
    group_split(X1) |>
    future_map_dfr(~mutate(., gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
			   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) |>
		   select(chr = X1, start = X4, end = X5, gene_id, gene_name))


# Standard ASE
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

write_tsv(ase_df, "./ase_data_raw.tsv")

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
    group_split(stim) |>
    future_map_dfr(~mutate(., p_value = map2_dbl(ref_count, ref_count + alt_count, 
						~binom.test(.x, .y, p = .5, alternative = "two.sided")$p.value)))

ase_res <- ase_res |> 
    mutate(q_value = qvalue(p_value)$qvalues,
	   p_bonferroni = p.adjust(p_value, method = "bonferroni"))

write_tsv(ase_res, "./ase_data.tsv")

plan(sequential)


# Dynamic ASE
ase_data_dyn <- 
    "./results_glm/glm_res_df.tsv" |>
    read_tsv(col_types = "cfccddd") 

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
ldsc_genes <- 
    ase_dyn_annotated |>
    filter(grepl("^Day 0-", stim), 
	   !is.na(gene_id)) |>
    separate_rows(gene_id, gene_name, sep = "/") |>
    mutate(fdr = p.adjust(p_dev, method = "fdr")) |>
    group_by(stim, gene_id, gene_name) |>
    summarise(is_signif = any(fdr <= 0.1)) |>
    ungroup()

write_tsv(ldsc_genes, "./ldsc_genes.tsv")

ldsc_genes_hg38 <- ldsc_genes |>
    left_join(gene_annot, join_by(gene_id, gene_name)) |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$"))

write_tsv(ldsc_genes_hg38, "./ldsc_genes_hg38.tsv")

# Use liftOver to convert annotations to hg19
dir.create("temp")
chain <- "/reference_databases/ReferenceGenome/liftover_chain/hg38/hg38ToHg19.over.chain.gz"
bedfile_38 <- "./temp/ldsc_genes_hg38.bed"
bedfile_19 <- "./temp/ldsc_genes_hg19.bed"
fail <- "./temp/fail.txt"

bed38 <- ldsc_genes_hg38 |>
    distinct(chr, start, end, gene_id) |>
    mutate(chr = factor(chr, levels = paste0("chr", c(1:22, "X"))),
	   start = start - 1L) |>
    arrange(chr, start, end)

write_tsv(bed38, bedfile_38, col_names = FALSE)

command <- glue("liftOver {bedfile_38} {chain} {bedfile_19} {fail}")
system(command)

bed19 <-
    read_tsv(bedfile_19, col_names = c("chr", "start", "end", "gene_id")) |>
    mutate(start = start + 1L)

ldsc_genes_hg19 <-
    ldsc_genes |>
    mutate(gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    left_join(bed19, join_by(gene_id))

write_tsv(ldsc_genes_hg19, "./ldsc_genes_hg19.tsv")
