library(tidyverse)
library(readxl)
library(rvest)
#library(furrr)

# Functions to read data
parse_pvalues <- function(p) {

    p %>%
	gsub("\\s", "", .) %>%
	sub("^([0-9.]+)[^0-9.]10[^0-9.](\\d+)$", "\\1e-\\2", .) %>%
	parse_number()
}

read_langefeld <- function(tab_path, sheet_number) {
    read_excel(tab_path, sheet_number, skip = 2) %>%
    select(snp_id = 1, 
	   chr = 2, 
	   pos = 3, 
	   gene_region = 4, 
	   region_rank = 5, 
	   ref_allele = 6, 
	   p = `P-value`,
	   or = `OR (95% CI)`,
	   p_stepwise = `Regional Stepwise  P-value`) %>% 
    filter(!is.na(pos)) %>%
    mutate(snp_id = gsub("\\s|[a-zA-Z,]+$", "", snp_id),  
	   p = parse_pvalues(p),
	   p_stepwise = parse_pvalues(p_stepwise),
	   or = parse_number(or)) %>%
    select(snp_id, chr, pos, gene_region, region_rank, ref_allele, p, p_stepwise, or)
}

read_bentham <- function(url_address) {

    read_html(url_address) %>%
    html_node("table") %>%
    html_table(header = TRUE, fill = TRUE) %>%
    select(snp_id = 1, chr = 2, pos = 3, locus = 4, p = 11, or = 12) %>%
    slice(-1) %>%
    mutate(snp_id = gsub("\\s|[a-zA-Z,]+$", "", snp_id),
	   pos = parse_number(pos),
	   p = parse_pvalues(p),
	   or = parse_number(or))
}

#get_var_info <- function(dbSNP_id) {
#    
#    "https://www.ncbi.nlm.nih.gov/snp/%s" %>%
#    sprintf(dbSNP_id) %>%
#    read_html() %>%    
#    html_nodes("div") %>%
#    html_nodes(".summary-box .usa-width-one-half") %>%
#    html_nodes("dd") %>%
#    html_text() %>%
#    .[4] %>%
#    trimws() %>%
#    gsub("^([^\n]+).*$", "\\1", .)
#}

# Import Langefeld data
tier1_ea <- read_langefeld("./sle_variants/langefeld_tableS2.xlxs", 1)
tier2_ea <- read_langefeld("./sle_variants/langefeld_tableS2.xlxs", 2)
tier3_ea <- read_langefeld("./sle_variants/langefeld_tableS2.xlxs", 3)

langefeld_df <- bind_rows("1" = tier1_ea, "2" = tier2_ea, "3" = tier3_ea, .id = "tier") %>%
    arrange(tier, region_rank, p)

langefeld_df %>% filter(grepl("imm", snp_id))

# Process Bentham et al. data
bentham_df <- read_bentham("https://www.nature.com/articles/ng.3434/tables/1")

bentham_supp <- 
    tribble(~snp_id, ~chr, ~pos, ~locus, ~p, ~or,
	    "rs10753074", 1, 173346343, "TNFSF4", 5.82e-12, 1.21,
	    "rs6736175", 2, 191946322, "STAT4", 9.17e-17, 1.24,
	    "rs9273076", 6, 32612301, "HLA-DQA1", 7.54e-13, 1.30,
	    "rs114092478", 6, 32682135, "HLA-DQB1", 2.9e-93, 2.0,
	    "rs114090659", 6, 30940989, "DPCR", 5.81e-92, 1.99,
	    "rs74290525", 6, 31835162, "SLC44A4", 1.12e-12, 2.06,
	    "rs3757387", 7, 128576086, "IRF5", 1.14e-48, 1.45)

# TLR7 variant (see review by Teruel & Arlacon-Riquelme)
deng_df <- 
    tibble(snp_id = "rs3853839",
	   chr = "X", 
	   pos = 12889539, 
	   locus = "TLR7", 
	   p = 2e-9,
	   or = 1.12)

# Combine studies
langefeld_info <- langefeld_df %>%
    filter(is.na(p_stepwise) | p_stepwise < 0.001) %>%
    select(chr, snp_id, pos, or) %>%
    mutate(chr = sub("^(\\d+).*$", "chr\\1", chr))

bentham_info <- bentham_df %>%
    select(chr, snp_id, pos, or) %>%
    mutate(chr = paste0("chr", chr)) %>%
    anti_join(langefeld_info, by = c("chr", "snp_id", "pos"))

bentham_supp_info <- bentham_supp %>%
    select(chr, snp_id, pos, or) %>%
    mutate(chr = paste0("chr", chr)) %>%
    anti_join(langefeld_info, by = c("chr", "snp_id", "pos"))

deng_info <- deng_df %>%
    select(chr, snp_id, pos, or)

info_df <- 
    bind_rows(langefeld = langefeld_info,
	      bentham = bentham_info,
	      bentham = bentham_supp_info,
	      deng = deng_info,
	      .id = "study") %>%
    mutate(chr = factor(chr, levels = paste0("chr", c(1:22, "X")))) %>%
    arrange(chr, pos)

#info_out <- info_df %>%
#    mutate(type = future_map_chr(snp_id, get_var_info))

write_tsv(info_df, "./sle_variants/sle_variants.tsv")
