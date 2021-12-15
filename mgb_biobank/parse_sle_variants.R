library(tidyverse)
library(readxl)
library(rvest)

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

# Import Langefeld data
tier1_ea <- read_langefeld("./sle_variants/langefeld_tableS2.xlxs", 1)
tier2_ea <- read_langefeld("./sle_variants/langefeld_tableS2.xlxs", 2)
tier3_ea <- read_langefeld("./sle_variants/langefeld_tableS2.xlxs", 3)

langefeld_tab <- bind_rows("1" = tier1_ea, "2" = tier2_ea, "3" = tier3_ea, .id = "tier") %>%
    arrange(tier, region_rank, p)

# Process Bentham et al. data
bentham_df <- read_bentham("https://www.nature.com/articles/ng.3434/tables/1")


# Combine Langefeld and Bentham
langefeld_bed <- langefeld_tab %>%
    filter(is.na(p_stepwise) | p_stepwise < 0.001) %>%
    select(snp_id, chr, pos) %>%
    mutate(chr = sub("^(\\d+).*$", "chr\\1", chr),
	   chr = factor(chr, levels = paste0("chr", c(1:22, "X"))),
	   end = pos + 1L) %>%
    select(chr, start = pos, end, snp_id) %>%
    arrange(chr, start) 

bentham_bed <- bentham_df %>%
    mutate(chr = paste0("chr", chr),
	   chr = factor(chr, levels = paste0("chr", c(1:22, "X"))),
	   end = pos + 1L) %>%
    select(chr, start = pos, end, snp_id)
   
final_bed <- 
    bind_rows("langefeld" = langefeld_bed, "bentham" = bentham_bed, .id = "study") %>%
    distinct(chr, start, end, snp_id, .keep_all = TRUE) %>%
    unite("id", c(snp_id, study), sep = "-") %>%
    select(chr, start, end, id) %>%
    arrange(chr, start)

write_tsv(final_bed, "./sle_variants/sle_langefeld_bentham.bed", col_names = FALSE)

langefeld_info <- langefeld_tab %>%
    select(chr, snp_id, pos, or) %>%
    mutate(chr = sub("^(\\d+).*$", "chr\\1", chr),
	   chr = factor(chr, levels = paste0("chr", c(1:22, "X"))))

bentham_info <- bentham_df %>%
    select(chr, snp_id, pos, or) %>%
    mutate(chr = paste0("chr", chr),
	   chr = factor(chr, levels = paste0("chr", c(1:22, "X"))))

info_df <- 
    bind_rows(langefeld = langefeld_info,
	      bentham = bentham_info, 
	      .id = "study") %>%
    arrange(chr, pos) %>%
    select(-pos)
   
write_tsv(info_df, "./sle_variants/sle_langefeld_bentham_info.tsv")
