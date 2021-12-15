library(tidyverse)

vars_df <- 
    "./sle_variants/sle_langefeld_bentham_hg38.bed" %>%
    read_tsv(col_names = c("chr", "start", "end", "info")) %>% 
    separate(info, c("var_id", "study"), sep = "-")

langefeld_df <- filter(vars_df, study == "langefeld") %>%
    select(chr, var_id, pos = start, study)

bentham_df <- filter(vars_df, study == "bentham") %>%
    select(chr, var_id, pos = start, study)

pairs_df <- right_join(langefeld_df, bentham_df, by = "chr") %>%
    filter(!is.na(var_id.x))

vcf <- "./sle_variants/sle.MGB.vcf" %>%
    read_tsv(comment = "##") %>%
    select(chr = 1, pos = 2, matches("^[0-9]")) %>%
    inner_join(select(vars_df, chr, pos = start)) %>%
    pivot_longer(-(1:2), names_to = "sample_id", values_to = "genotype") %>%
    group_by(chr, pos) %>%
    filter(!(any(genotype == "./."))) %>%
    ungroup()

vcf_dose <- vcf %>%
    mutate(dose = case_when(genotype == "0|0" ~ 0L,
			    genotype == "0|1" ~ 1L,
			    genotype == "0|2" ~ 1L,
			    genotype == "1|0" ~ 1L,
			    genotype == "1|1" ~ 2L)) %>%
    select(-genotype)

cor_df <- pairs_df %>%
    left_join(vcf_dose, by = c("chr", "pos.x" = "pos")) %>%
    filter(!is.na(sample_id)) %>%
    left_join(vcf_dose, by = c("chr", "sample_id", "pos.y" = "pos")) %>%
    group_by(chr, study.x, var_id.x, study.y, var_id.y) %>%
    summarise(r2 = cor(dose.x, dose.y)^2) %>%
    ungroup()

out <- cor_df %>%
    mutate(chr = factor(chr, levels = paste0("chr", c(1:22, "X")))) %>%
    left_join(select(langefeld_df, var_id, pos), by = c("var_id.x" = "var_id")) %>% 
    arrange(chr, pos) %>%
    mutate(var_id.y = fct_inorder(var_id.y)) %>%
    mutate(region = group_indices(., var_id.y)) %>%
    select(chr, region, langefeld = var_id.x, bentham = var_id.y, pos, r2)
   
write_tsv(out, "./sle_variants/sle_ld.tsv")
