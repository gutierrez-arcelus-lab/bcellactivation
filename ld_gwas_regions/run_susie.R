library(tidyverse)
library(susieR)
library(glue)
library(ggrepel)

risk_var <- "rs7582694"

ld_vcf <- 
    "./data/chr2:190970120-192970120.vcf.gz" |>
    data.table::fread(skip = "#CHROM") |>
    as_tibble() |>
    unite("ID", c(ID, REF, ALT), sep = "-")

ld_plink <- 
    "./data/chr2:190970120-192970120_r.ld" |>
    data.table::fread() |>
    data.matrix()

rownames(ld_plink) <- colnames(ld_plink) <- ld_vcf$ID

sample_size <- 4036 + 6959

langefeld <- 
    "./data/langefeld_stat4.tsv" |>
    data.table::fread() |>
    mutate(chrom = "2") |>
    select(chrom, pos, rsid, other_allele = alleleA, effect_allele = alleleB, p, beta, se) |>
    filter(between(pos, pos - 0.75e5, pos + 0.75e5)) |>
    as_tibble() |>
    mutate(snp_id = glue("{rsid}-{other_allele}-{effect_allele}")) |>
    filter(snp_id %in% colnames(ld_plink))

ld_plink_final <- ld_plink[langefeld$snp_id, langefeld$snp_id]

langefeld_final <- 
    langefeld |>
    mutate(snp_id = factor(snp_id, levels = colnames(ld_plink))) |>
    arrange(snp_id) |>
    select(chrom, pos, snp_id, p, beta, se) |>
    mutate(z = beta/se)

condz <- kriging_rss(langefeld_final$z, ld_plink_final, n = sample_size, r_tol = 1e-04)
png("./plots/susie_kriging_rss_STAT4.png")
condz$plot
dev.off()

condz_df <- 
    condz[[2]] |>
    as_tibble() |>
    bind_cols(select(langefeld_final, snp_id))

#condz_plot <- 
#    ggplot(condz_df, aes(x = condmean, y = z)) +
#    geom_abline() +
#    geom_point(size = .5) +
#    geom_text_repel(data = top_n(condz_df, 3, abs(z_std_diff)),
#		    aes(x = condmean, y = z, label = snp_id),
#		    min.segment.length = 0.1) +
#    theme_bw() +
#    labs(x = "Expected value", y = "Observed value")
#
#ggsave("./plots/susie_kriging_rss_STAT4.png", condz_plot, width = 5, height = 5)

# if it does not converge, remove problematic SNPs and repeat until convergence
iter <- 0L
fit <- list()
fit$converged <- FALSE
summ_stats <- langefeld_final
ldmat <- ld_plink_final

while ( !fit$converged & iter <= 30 ) {

    if ( iter > 0 ) { 
    
	condz <- kriging_rss(summ_stats$z, 
			     ldmat, 
			     n = sample_size, 
			     r_tol = 1e-04)

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

coverage_df <- 
    tibble(coverage = scales::percent(round(fit$sets$coverage, 3)),
	   cs = paste0("L", fit$sets$cs_index))

cs_df <- 
    enframe(fit$sets$cs, "cs", "rowid") |>
    unnest(cols = rowid) |>
    left_join(coverage_df, by = "cs")

pip_df <- 
    enframe(fit$pip, "ID", "pip") |>
    separate(ID, c("rsid", "ref", "alt"), sep = "-", convert = TRUE) |>
    rowid_to_column() |>
    left_join(cs_df, join_by(rowid)) |>
    mutate(cs = ifelse(!is.na(cs), paste0(cs, " (", coverage, ")"), NA))

write_tsv(pip_df, "./data/susie_stat4.tsv")
