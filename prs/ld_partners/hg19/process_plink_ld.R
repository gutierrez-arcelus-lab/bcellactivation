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
	mutate(idmin = pmin(id1, id2),
	       idmax = pmax(id1, id2)) |>
	distinct(idmin, idmax, r2) |>
	filter(idmin != idmax, r2 >= .8) |>
	filter(idmin %in% snp_info$rsid & idmax %in% snp_info$rsid)

    if ( nrow(ldmat_filt) == 0 ) {
	return(tibble(sentinel = s_snp, snp2 = NA, r2 = NA)) 
    }

    out <- 
	ldmat_filt |>
	select(snp1 = idmin, snp2 = idmax, r2) |>
	filter(snp1 == s_snp | snp2 == s_snp) |>
	mutate(sentinel = ifelse(snp1 == s_snp, snp1, snp2),
	       snp2 = ifelse(snp2 == s_snp, snp1, snp2)) |>
	select(sentinel, snp2, r2) |>
	arrange(desc(r2))

    if ( nrow(out) == 0 ) {
	return(tibble(sentinel = s_snp, snp2 = NA, r2 = NA)) 
    } else {
	return(out)
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

write_tsv(res_df, "./data/Khunsriraksakul_sentinel_r2.tsv")
