library(tidyverse)

s_snp <- commandArgs(TRUE)[1]

out_dir <- file.path(system("echo $TEMP_WORK", intern = T), "vcf/prs/hg19/pt_5e-4")

vcf <- 
    file.path(out_dir, "%s.vcf.gz") |>
    sprintf(s_snp) |>
    read_tsv(comment = "##")  |>
    unite("ID", c(ID, REF, ALT), sep = "_") |>
    select(chr = `#CHROM`, POS, ID)

plink <- 
    file.path(out_dir, "%s.ld") |>
    sprintf(s_snp) |>
    read_delim(delim = " ", col_names = FALSE)

if ( (ncol(plink) %% nrow(plink) == 1) && all(is.na(plink[, ncol(plink)])) ) {

    plink <- plink[, -ncol(plink)]
}

colnames(plink) <- vcf$ID
plink <- add_column(plink, snp1 = vcf$ID, .before = 1)
valid_snps <- names(which(sapply(plink, function(x) !all(is.na(x)))))

out <- plink |>
    filter(snp1 %in% valid_snps) |>
    select(snp1, all_of(valid_snps)) |>
    pivot_longer(-snp1, names_to = "snp2", values_to = "r2") |>
    filter(grepl(sprintf("^%s_", s_snp), snp1), 
	   r2 >= 0.8) |>
    arrange(desc(r2)) |>
    mutate(var_type = case_when(snp1 == snp2 ~ "lead",
				TRUE ~ "tag")) |>
    select(snp2, var_type, r2) |>
    left_join(vcf, join_by(snp2 == ID)) |>
    separate(snp2, c("rsid", "ref", "alt"), sep = "_") |>
    select(chr, pos = POS, rsid, var_type, ref, alt, r2)

out_file <- file.path(out_dir, s_snp) |>
    paste0("_tagvariants.tsv")

write_tsv(out, out_file)
