library(edgeR)
library(splines)
library(tidyverse)

# Functions

## filter and normalize to TMM
filt_norm <- function(y) {

    keep_y <- rowMeans(cpm(y) > 1) >= 0.7 
    y <- y[keep_y, , keep.lib.sizes = FALSE]

    y <- calcNormFactors(y)
    y <- estimateDisp(y)
    
    return(y)
}


# data
expr_df <- read_rds("./data/expression_transcripts.rds") %>%
    mutate(gene_id = sub("\\.\\d+$", "", gene_id))

counts <- expr_df %>%
    select(-tpm, -gene_name, -gene_id) %>%
    pivot_wider(names_from = id, values_from = count) %>%
    mutate(tx_id = sub("\\.\\d+$", "", tx_id)) %>%
    column_to_rownames("tx_id") %>%
    data.matrix()

sample_decode <- read_tsv("./data/sample_decode.tsv") %>%
    separate(sample_name, c("name", "stim", "time"), sep = "_", remove = FALSE) %>%
    mutate(time = factor(time, levels = str_sort(unique(time), numeric = TRUE))) %>%
    arrange(stim, time, name)

counts_sorted <- counts[, sample_decode$sample_id]

pooled <- sumTechReps(counts_sorted, ID=sample_decode$sample_name)

# Use Unstim 0hr as 0hr for all stims
matrix_0hrs <- pooled[ , grep("_0hrs$", colnames(pooled), value = TRUE)] %>%
    as_tibble(rownames = "gene_id") %>%
    pivot_longer(-gene_id, names_to = "sample_name") %>%
    separate(sample_name, c("id", "stim", "time"), sep = "_") %>%
    select(-stim) %>%
    nest(data = everything()) %>%
    expand_grid(stim = unique(sample_decode$stim)) %>%
    filter(stim != "Unstim") %>%
    unnest(cols = c(data)) %>%
    unite("sample_name", c(id, stim, time), sep = "_") %>%
    pivot_wider(names_from = sample_name, values_from = value) %>%
    column_to_rownames("gene_id") %>%
    data.matrix() %>%
    .[rownames(pooled), ]



counts_full_tmp <- cbind(pooled, matrix_0hrs) %>%
    as_tibble(rownames = "tx_id") %>%
    pivot_longer(-tx_id, names_to = "sample_name")

# break the pipe chain here
# because separate() on the whole dataset was taking too long

sample_df <- distinct(counts_full_tmp, sample_name) %>%
    separate(sample_name, c("id", "stim", "time"), sep = "_", remove = FALSE) %>%
    mutate(time = factor(time, levels = str_sort(unique(time), numeric = TRUE)))

counts_full <- counts_full_tmp %>%
    left_join(sample_df, by = "sample_name") %>%
    arrange(stim, time, id) %>%
    select(tx_id, sample_name, value) %>%
    pivot_wider(names_from = "sample_name", values_from = "value") %>%
    column_to_rownames("tx_id") %>%
    data.matrix()

stim_matrices <- sample_decode %>%
    distinct(stim) %>%
    filter(! stim %in% c("Unstim", "IL4")) %>%
    mutate(mat = map(stim, ~counts_full[, grepl(sprintf("_%s_", .x), colnames(counts_full))]),
	   tp = map(mat, ~sub("^[^_]+_[^_]+_(\\d+hrs)$", "\\1", colnames(.x))),
	   dge = map2(mat, tp, ~DGEList(counts = .x, group = .y)),
	   dgefilt = map(dge, filt_norm),
	   design = map(tp, function(x) {
			    hours = parse_number(x)
			    X = ns(hours, df = 3)
			    model.matrix(~X)}),
	   fit = map2(dgefilt, design, ~glmQLFit(.x, .y, robust = TRUE)),
	   fit = map(fit, ~glmQLFTest(.x, coef = 2:4)),
	   logcpm_obs = map2(dgefilt, fit, ~cpm(.x, log = TRUE, prior.count = .y$prior.count)),
	   logcpm_fit = map(fit, ~cpm(.x, log = TRUE)))

obs_cpm_df <- stim_matrices %>%
    select(stim, logcpm_obs) %>%
    mutate(logcpm_obs = map(logcpm_obs, ~as_tibble(.x, rownames = "tx_id") %>%
			    pivot_longer(-tx_id, names_to = "sample_name", values_to = "cpm"))) %>%
    unnest(cols = c(logcpm_obs))

fit_cpm_df <- stim_matrices %>%
    select(stim, logcpm_obs, logcpm_fit) %>%
    mutate(logcpm_fit = map2(logcpm_fit, logcpm_obs, ~`colnames<-`(.x, colnames(.y)))) %>%
    select(-logcpm_obs) %>%
    mutate(logcpm_fit = map(logcpm_fit, ~as_tibble(.x, rownames = "tx_id") %>%
			    pivot_longer(-tx_id, names_to = "sample_name", values_to = "fit"))) %>%
    unnest(cols = c(logcpm_fit))

genes_df <- distinct(expr_df, gene_id, gene_name, tx_id) %>%
    mutate(tx_id = sub("\\.\\d+$", "", tx_id))

cpm_df <- left_join(obs_cpm_df, fit_cpm_df) %>%
    left_join(genes_df, by = "tx_id") %>%
    select(-stim) %>%
    left_join(sample_df, by = "sample_name") %>%
    select(id, stim, time, gene_id, gene_name, tx_id, cpm, fit)

results_df <- stim_matrices %>%
    select(stim, fit) %>%
    mutate(toptag = map(fit, ~topTags(.x, n = Inf) %>% 
			as.data.frame() %>%
			as_tibble(rownames = "tx_id") %>%
			left_join(genes_df, by = "tx_id") %>%
			select(gene_id, gene_name, tx_id, everything()))) %>%
    select(stim, toptag) %>%
    unnest(cols = c(toptag))

write_tsv(cpm_df, "./data/edger_transcript_cpm_fit.tsv")
write_tsv(results_df, "./data/edger_de_transcripts.tsv")
