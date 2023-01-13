library(coloc)
library(GenomicRanges)
library(seqminer)
library(tidyverse)
library(furrr)

# functions
import_eqtl <- function(ftp_path, region, header, genes) {
 
    fetch <- safely(function(f, r) {
	tabix.read.table(tabixFile = f, 
			 tabixRange = r) |>
	as_tibble()
    })

    fetch_table <- fetch(f = ftp_path, r = region) 
    error <- fetch_table$error
	
    i <- 1
    while ( !is.null(error) & i < 30 ) {
	Sys.sleep(10)
	fetch_table <- fetch(f = ftp_path, r = region) 
	error <- fetch_table$error
	cat("Retrying...\n")
	i <- i + 1
    }
    
    if ( !is.null(fetch_table$error) ) stop("Can't retrieve dataset")

    out <- as_tibble(fetch_table$result)
    colnames(out) <- header

    eqtl <- filter(out, gene_id %in% genes) |>
	mutate(rsid = sub("\\\r$", "", rsid)) |>
	filter(rsid != "NA") |>
	add_count(molecular_trait_id, gene_id, rsid) |>
	filter(n == 1) |>
	select(-n) |>
	add_count(molecular_trait_id, gene_id, type, position) |>
	filter(n == 1) |>
	select(-n)
    
    eqtl_cat_df <- eqtl |>
	select(gene_id, molecular_trait_id, rsid, position, pvalue) 
    
    eqtl_cat_df
}

# Directoy
rundir <- "/lab-share/IM-Gutierrez-e2/Public/vitor/lupus/colocalization"

procs <- 8

# eQTL Catalogue column names
column_names <- read_lines(file.path(rundir, "eqtl_column_names.txt"))

# eQTL catalogue URLs
ftps <- file.path(rundir, "data/coloc_input/ftp_urls.txt") |>
    read_lines()

study_ids <- basename(ftps) |> 
    sub("\\.all\\.tsv\\.gz", "", x = _)

names(ftps) <- study_ids

# regions
region_df <- file.path(rundir, "data/coloc_input/region.tsv") |>
    read_tsv() |>
    filter(region %in% c(5, 22)) |>
    add_column(study = "bentham", .before = 1)

region_df_2 <- file.path(rundir, "data/coloc_input/region_langefeld.tsv") |>
    read_tsv() |>
    filter(region %in% c(22, 31)) |>
    add_column(study = "langefeld", .before = 1)

region_all_df <- bind_rows(region_df, region_df_2) |>
    extract(coord, c("chr", "start", "end"), "(\\d):(\\d+)-(\\d+)", convert = TRUE) |>
    group_by(chr) |>
    summarise(start = min(start),
	      end = max(end)) |>
    ungroup() |>
    mutate(coord = paste0(chr, ":", start, "-", end)) |>
    select(coord)

# tissues of relevance
ftp_df <- "./data/coloc_input/eqtl_catalogue_paths.tsv" |>
    read_tsv() |>
    mutate(study_label = basename(ftp_path),
	   study_label =  sub("\\.all\\.tsv\\.gz", "", x = study_label)) |>
    select(study = study_label, author = study, tissue = tissue_label,
	   condition = condition_label, qtl_group, method = quant_method)

tissues <- c("macrophage", "monocyte", "neutrophil", "CD4+ T cell", "CD8+ T cell",
	     "B cell", "LCL", "T cell", "blood", "Tfh cell", "Th1 cell", "Th2 cell",
	     "Treg naive", "Treg memory", "CD16+ monocyte", "NK cell")

studies <- filter(ftp_df, tissue %in% tissues) |>
    pull(study)

design_df <- tibble(eqtl_study = study_ids, ftp = ftps) |>
    filter(eqtl_study %in% studies) |>
    cross_join(region_all_df)

genes_df <- read_tsv("./data/coloc_input/gene_annots_langefeld.tsv") |>
    filter(gene_name %in% c("IKBKE", "IKZF1"))

# Run analysis
plan(multisession, workers = procs)

coloc_res <- design_df |>
    mutate(dat = future_map2(ftp, coord, 
			     safely(function(x, y) import_eqtl(ftp_path = x, 
							       region = y, 
							       header = column_names,
							       genes = genes_df$gene_id))))

coloc_res |>
    mutate(error = map(dat, 2)) |>
    pull(error) |>
    keep(function(x) !is.null(x))

coloc_out <- coloc_res |>
    mutate(dat = map(dat, 1)) |>
    select(eqtl_study, dat) |>
    unnest(cols = dat)

out_name <- file.path(rundir, "results", "eqtlcatalog_ikbke_ikzf1.tsv") 

write_tsv(coloc_out, out_name)
