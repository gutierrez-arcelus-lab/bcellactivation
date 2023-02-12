library(edgeR)
library(splines)
library(annotables)
library(tidyverse)


# Functions
run_edger <- function(stim_i) {

    counts_stim <- expr_df |>
	filter(stim %in% c("unstday0", stim_i)) |>
	unite("sample_id", c(sample_id, stim), sep = "_") |>
	select(-tpm, -gene_name) |>
	pivot_wider(names_from = sample_id, values_from = counts) |>
	column_to_rownames("gene_id")

    reps_df <- tibble(sample_id = names(counts_stim)) |>
	mutate(final_id = sub("^([^.]+)\\.[123]_(.+)$", "\\1_\\2", sample_id))

    pooled <- sumTechReps(counts_stim, ID = reps_df$final_id)
    groups_stim <- sub("^\\d+_(.+)$", "\\1", colnames(pooled))
    colnames(pooled) <- sub("^(\\d+)_.+$", "\\1", colnames(pooled))

    subject <- factor(colnames(pooled))
    stims <- factor(groups_stim, levels = c("unstday0", stim_i))
    design <- model.matrix(~subject+stims)

    y <- DGEList(counts = pooled, group = stims)
    keep_stim <- filterByExpr(y)
    y <- y[keep_stim, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit)

    cpm_obs <- cpm(y, log = FALSE, prior.count = qlf$prior.count)
    colnames(cpm_obs) <- paste(colnames(cpm_obs), groups_stim, sep = "_")

    cpm_df <- cpm_obs |>
	as.data.frame() |>
	rownames_to_column("gene_id") |>
	pivot_longer(-gene_id, names_to = "sample_id", values_to = "cpm") |>
	separate(sample_id, c("sample_id", "stim"), sep = "_") |>
	mutate(stim = factor(stim, levels = levels(stims)))

    res <- topTags(qlf, n = nrow(cpm_obs)) |>
	as.data.frame() |>
	as_tibble(rownames = "gene_id")

    list(cpm.data = cpm_df, results = res, fit = qlf)
}


# data
expr_df <- read_tsv("./quantifications_genes.tsv") |>
    mutate(gene_id = sub("\\.\\d+$", "", gene_id),
	   sample_id = as.character(sample_id))

gene_ids <- distinct(expr_df, gene_id, gene_name)

# Run analysis
res_df <- tibble(stim_test = c("BCR", "TLR7", "DN2")) |>
    mutate(dat = map(stim_test, run_edger))

cpm_final <- res_df |>
    mutate(cpm = map(dat, "cpm.data")) |>
    select(-dat) |>
    unnest(cols = cpm) |>
    left_join(gene_ids) |>
    select(stim_test, gene_id, gene_name, everything())

results_final <- res_df |>
    mutate(results = map(dat, "results")) |>
    select(-dat) |>
    unnest(cols = results) |>
    left_join(gene_ids) |>
    select(stim_test, gene_id, gene_name, everything())

write_tsv(cpm_final, "./edger_cpm.tsv")
write_tsv(results_final, "./edger_results.tsv")

# GO

fit_df <- res_df |>
    mutate(fit = map(dat, "fit")) |>
    select(-dat)

entrez <- grch38 |> 
    group_by(ensgene) |> 
    slice(1) |>
    ungroup() |>
    select(ensgene, entrez, symbol)

run_go <- function(fitobj) {
 
    genes <- rownames(fitobj)
    genes <- genes[genes %in% entrez$ensgene]

    fit <- fitobj[genes, ]

    genes_dat <- inner_join(tibble(ensgene = rownames(fit)), entrez)

    go <- goana(fit, geneid = genes_dat$entrez, species = "Hs")

    topGO(go, sort = "up", ontology = "BP", number = Inf) |> 
	as_tibble(rownames = "goid") |>
	filter(! Term %in% c("biological_process", "cellular process", 
			     "metabolic process", "cellular metabolic process", 
			     "biological regulation"))
}


go_res <- fit_df |>
    mutate(go = map(fit, run_go)) |>
    select(stim = stim_test, go) |>
    unnest(cols = go)

top_diff_paths <- go_res |>
    group_by(goid, Term) |>
    filter(all(c("BCR", "TLR7", "DN2") %in% stim)) |>
    filter(any(P.Up < 0.01)) |>
    ungroup() |>
    group_split(stim) |>
    map_df(~rowid_to_column(., "i")) |>
    select(stim, i, goid, Term) |>
    pivot_wider(names_from = stim, values_from = i) |>
    mutate(d1 = abs(BCR - DN2), 
	   d2 = abs(BCR - TLR7),
	   d3 = abs(DN2 - TLR7)) |>
    select(goid, Term, d1, d2, d3) |>
    pivot_longer(d1:d3, names_to = "d", values_to = "value") |>
    group_by(goid, Term) |>
    slice(which.max(value)) |>
    ungroup() |>
    arrange(desc(value)) |>
    slice(1:24) |>
    select(goid, Term)

inner_join(top_diff_paths, go_res, multiple = "all") |>
    write_tsv("./top_diff_go.tsv")


