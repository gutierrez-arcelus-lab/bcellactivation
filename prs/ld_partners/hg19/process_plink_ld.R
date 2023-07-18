library(tidyverse)
library(furrr)

out_dir <- file.path(system("echo $TEMP_WORK", intern = T), "vcf/prs/hg19")

snp_info <-
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Khunsriraksakul/Khunsriraksakul_summstats_mamt.txt" |>
    read_tsv() |>
    janitor::clean_names() |>
    select(chr, pos = pos_hg19, rsid, ref = effect_allele, alt = other_allele)

sentinels <- read_tsv("./data/regions.tsv", col_names = c("region", "rsid"))

main <- function(s_snp) {
    
    vcf <- 
	file.path(out_dir, "%s.filt.vcf.gz") |>
	sprintf(s_snp) |>
	read_tsv(comment = "##") 

    plink <- 
	file.path(out_dir, "%s.ld") |>
	sprintf(s_snp) |>
	read_delim(delim = " ", col_names = FALSE)

    if ( (ncol(plink) %% nrow(plink) == 1) && all(is.na(plink[, ncol(plink)])) ) {

	plink <- plink[, -ncol(plink)]
    }

    valid_snps <- which(sapply(plink, function(x) !all(is.na(x))))

    plink_filt <- plink[valid_snps, valid_snps]
    ldmat <- data.matrix(plink_filt)
    vcf <- slice(vcf, valid_snps)
    rownames(ldmat) <- vcf$ID
    colnames(ldmat) <- vcf$ID

    ldmat_filt <- 
	ldmat |>
	as_tibble(rownames = "id1") |>
	pivot_longer(-id1, names_to = "id2", values_to = "r2") |>
	filter(id1 == s_snp, 
	       id1 != id2,
	       id2 %in% snp_info$rsid) |>
	arrange(desc(r2)) |>
	select(sentinel = id1, snp2 = id2, r2)

    if ( nrow(ldmat_filt) == 0 ) {
	return(tibble(sentinel = s_snp, snp2 = NA, r2 = NA)) 
    } else {
	return(ldmat_filt)
    }

}

# run analysis
plan(multisession, workers = availableCores())

res <- future_map(sentinels$rsid, safely(~main(s_snp = .x)))

errors <- map(res, "error")

errors |>
	map(~!is.null(.)) |>
	unlist() |>
	which()

res_df <- res |>
    map("result") |>
    bind_rows()

sentinels |>
    select(sentinel = rsid) |>
    left_join(filter(res_df, r2 >= .8)) |>
    write_tsv("./data/Khunsriraksakul_sentinel_r2.tsv")

read_tsv("./data/Khunsriraksakul_sentinel_r2.tsv") |>
    select(-r2) |>
    pivot_longer(sentinel:snp2) |>
    drop_na() |>
    distinct() |>
    pull(value) |>
    write_lines("./data/sentinels_and_LDpartners_rsid.txt")
