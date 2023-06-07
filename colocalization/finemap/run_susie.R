# Use 48Gb of RAM


library(tidyverse)
library(susieR)
library(furrr)

# Functions
run_susie <- function(locus_name) {

    region_id <- regions$region[regions$locus == locus_name] 

    vcf <- sprintf("./data/%s.adj.vcf.gz", region_id) |>
	read_tsv(comment = "##") 

    plink <- sprintf("./data/%s.ld", region_id) |>
	read_delim(delim = " ", col_names = FALSE)

    plink_filt <- plink[, sapply(plink, function(x) !all(is.na(x)))]
    ldmat <- data.matrix(plink_filt)
    rownames(ldmat) <- vcf$ID
    colnames(ldmat) <- vcf$ID

    summ_stats <- summ_stats_all |>
	filter(locus == locus_name) |>
	unite("ID", c(chr, bp, ref, alt), sep = ":") |>
	filter(ID %in% vcf$ID) |>
	mutate(z = beta/se)

    #condz <- kriging_rss(summ_stats$z, ldmat, n = sample_size, r_tol = 1e-04)
    #png(sprintf("./plots/diagnostic_%s.png", locus_name))
    #condz$plot
    #dev.off()

    # if it does not converge, remove problematic SNPs and repeat until convergence
    iter <- 0L
    fit <- list()
    fit$converged <- FALSE
    while ( !fit$converged & iter <= 20 ) {

	if ( iter > 0 ) { 
	
	    condz <- kriging_rss(summ_stats$z, ldmat, n = sample_size, r_tol = 1e-04)

	    snp_rm <- as_tibble(condz$conditional_dist) |>
		rowid_to_column() |>
		arrange(desc(abs(z_std_diff))) |>
		slice(1) |>
		pull(rowid)

	    ldmat <- ldmat[-snp_rm, -snp_rm]
	    summ_stats <- slice(summ_stats, -snp_rm)
	}

	fit <- susie_rss(summ_stats$z, 
			 ldmat, 
			 n = sample_size, 
			 L = 10, 
			 coverage = 0.9, 
			 r_tol = 1e-05)

	iter <- iter + 1L
    }

    return(fit)
}

make_pip_df <- function(fit_obj) {

    if (! is.null(fit_obj$sets$coverage) ) {

	coverage_df <- 
	    tibble(coverage = scales::percent(round(fit_obj$sets$coverage, 3)),
		   cs = paste0("L", fit_obj$sets$cs_index))

	cs_df <- enframe(fit_obj$sets$cs, "cs", "rowid") |>
	    unnest(cols = rowid) |>
	    left_join(coverage_df, by = "cs")

	pip_df <- enframe(fit_obj$pip, "ID", "pip") |>
	    separate(ID, c("chr", "pos", "ref", "alt"), sep = ":", convert = TRUE) |>
	    rowid_to_column() |>
	    left_join(cs_df) |>
	    mutate(cs = ifelse(!is.na(cs), paste0(cs, " (", coverage, ")"), NA))

    } else {
	
	pip_df <- 
	    enframe(fit_obj$pip, "ID", "pip") |>
	    separate(ID, c("chr", "pos", "ref", "alt"), sep = ":", convert = TRUE) |>
	    mutate(cs = NA)
    }

    return(pip_df)
}

## Data
# Genomics windows
regions <- 
    "./data/regions_bentham_1mb.tsv" |>
    read_tsv(col_names = c("region", "locus")) |>
    extract(region, c("chr", "left", "right"), 
	    "(chr[^:]+):(\\d+)-(\\d+)", 
	    convert = TRUE, remove = FALSE)

# GWAS data
sample_size <- 4036 + 6959
Sys.setenv(VROOM_CONNECTION_SIZE = 500000L)
summ_stats_all <- read_tsv("./data/bentham_opengwas_1MbWindows_hg38_summstats.tsv")

# Analysis
plan(multisession, workers = availableCores())

susie_results <- regions$locus |>
    setNames(regions$locus) |>
    future_map(run_susie)

write_rds(susie_results, "./susie_results.rds")

pip_df <- map_dfr(susie_results, make_pip_df, .id = "locus")

write_tsv(pip_df, "./susie_results.tsv")
