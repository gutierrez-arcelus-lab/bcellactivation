library(tidyverse)

bed <- 
    "./sle_variants/sle_variants_hg38.bed" %>%
    read_tsv(col_names = c("chr", "start", "end")) %>% 
    select(chr, pos_hg38 = start)

vars_df <- read_tsv("./sle_variants/sle_variants.tsv") %>%
    add_column(pos_hg38 = bed$pos_hg38, .after = "pos")

langefeld_df <- filter(vars_df, study == "langefeld") %>%
    select(chr, snp_id, pos = pos_hg38, study)

bentham_df <- filter(vars_df, study == "bentham") %>%
    select(chr, snp_id, pos = pos_hg38, study)

pairs_df <- right_join(langefeld_df, bentham_df, by = "chr") %>%
    filter(!is.na(snp_id.x))

vcf <- "./sle_variants/sle.MGB.vcf" %>%
    read_tsv(comment = "##") %>%
    select(chr = 1, pos = 2, matches("^[0-9]")) %>%
    distinct(chr, pos, .keep_all = TRUE) %>%
    inner_join(select(vars_df, chr, pos = pos_hg38)) %>%
    pivot_longer(-(1:2), names_to = "sample_id", values_to = "genotype") %>%
    group_by(chr, pos) %>%
    filter(!(any(genotype == "./."))) %>%
    ungroup()

count(vcf, genotype)

vcf_dose <- vcf %>%
    mutate(dose = case_when(genotype == "0|0" ~ 0L,
			    genotype == "0|1" ~ 1L,
			    genotype == "1|0" ~ 1L,
			    genotype == "1|1" ~ 2L)) %>%
    select(-genotype)

cor_df <- pairs_df %>%
    left_join(vcf_dose, by = c("chr", "pos.x" = "pos")) %>%
    filter(!is.na(sample_id)) %>%
    left_join(vcf_dose, by = c("chr", "sample_id", "pos.y" = "pos")) %>%
    group_by(chr, study.x, snp_id.x, study.y, snp_id.y) %>%
    summarise(r2 = cor(dose.x, dose.y)^2) %>%
    ungroup()

ld_df <- cor_df %>%
    mutate(chr = factor(chr, levels = paste0("chr", c(1:22, "X")))) %>%
    left_join(select(langefeld_df, snp_id, pos), by = c("snp_id.x" = "snp_id")) %>% 
    arrange(chr, pos) %>%
    mutate(snp_id.y = fct_inorder(snp_id.y)) %>%
    group_by(snp_id.y) %>%
    mutate(region = cur_group_id()) %>%
    ungroup() %>%
    select(chr, region, langefeld = snp_id.x, bentham = snp_id.y, pos, r2)
   
write_tsv(ld_df, "./sle_variants/sle_ld.tsv")

ld_vars <- ld_df %>%
    filter(langefeld != bentham, r2 >= 0.6) %>%
    distinct(chr, bentham) %>%
    pivot_longer(-chr, names_to = "study", values_to = "snp_id") 

vars_indep <- vars_df %>%
    mutate(or = ifelse(or < 1L, 1L/or, or)) %>%
    anti_join(ld_vars, by = c("study", "chr", "snp_id")) %>%
    select(chr, pos = pos_hg38, or)

dosage_df <- vcf_dose %>%
    inner_join(vars_indep) %>%
    mutate(beta = log(or))

out <- dosage_df %>%
    group_by(sample_id) %>%
    summarise(het_score = sum(as.integer(dose == 1)),
	      het_score_wt = sum(as.integer(dose == 1) * beta))

write_tsv(out, "./sle_variants/scores.tsv")

