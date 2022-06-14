library(tidyverse)
library(tximport)
library(DRIMSeq)

# Cell type
cell_type <- commandArgs(TRUE)[1]

# Sample specs
samps <- file.path("./results", cell_type, "groups_file.txt") %>%
    read_delim(col_names = c("sample_id", "condition"), delim = " ") %>%
    mutate(condition = factor(condition))

# Salmon files
files <- file.path("/lab-share/IM-Gutierrez-e2/Public/vitor/scharer/salmon", samps$sample_id, "quant.sf")
names(files) <- samps$sample_id

# Import expression data with tximport
txi <- tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "scaledTPM")

# Remove EBV transcripts
txi$abundance <- txi$abundance[!grepl("^HHV4", rownames(txi$abundance)), ]
txi$counts <- txi$counts[!grepl("^HHV4", rownames(txi$counts)), ]
txi$length <- txi$length[!grepl("^HHV4", rownames(txi$length)), ]

cts <- txi$counts
cts <- cts[rowSums(cts) > 0, ]

range(colSums(cts)/1e6)

# Transcript to Gene table
tx_to_gene <- "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.txToGene.tsv" %>%
    read_tsv()

txdf <- select(tx_to_gene, TXNAME = transcript_id, GENEID = gene_id) %>%
    filter(TXNAME %in% rownames(cts)) %>%
    add_count(GENEID, name = "ntx")

# Create DRIMSeq data.frame
counts <- as.data.frame(cts) %>%
    rownames_to_column("TXNAME") %>%
    inner_join(txdf, .) %>%
    select(gene_id = GENEID, feature_id = TXNAME, starts_with("SRR"))

d <- dmDSdata(counts = as.data.frame(counts), samples = as.data.frame(samps))

# Filters
d_filt <- dmFilter(d, 
		   min_samps_feature_expr = 3, min_feature_expr = 10,
		   min_samps_feature_prop = 3, min_feature_prop = 0.1,
		   min_samps_gene_expr = 3, min_gene_expr = 10)

table(table(counts(d_filt)$gene_id))

# DTU analysis
design_full <- model.matrix(~condition, data = DRIMSeq::samples(d_filt))

set.seed(1)
d_filt <- dmPrecision(d_filt, design = design_full) %>%
    dmFit(design = design_full) %>%
    dmTest(coef = "conditionSLE")

res <- DRIMSeq::results(d_filt) %>%
    mutate_at(vars(pvalue, adj_pvalue), ~replace_na(., 1))

res_txp <- DRIMSeq::results(d_filt, level = "feature") %>%
    mutate_at(vars(pvalue, adj_pvalue), ~replace_na(., 1))

out_gene <- res %>%
    as_tibble() %>%
    filter(adj_pvalue < 0.05) %>%
    arrange(adj_pvalue) %>%
    left_join(distinct(tx_to_gene, gene_id, gene_name))

out_tx <- res_txp %>%
    as_tibble() %>%
    filter(adj_pvalue < 0.05) %>%
    arrange(adj_pvalue) %>%
    left_join(tx_to_gene, by = c("feature_id" = "transcript_id", "gene_id"))

out_file_gene <- file.path("./results", cell_type, "dtu_perGene.tsv") 
out_file_tx <- file.path("./results", cell_type, "dtu_perTx.tsv") 

write_tsv(out_gene, out_file_gene)
write_tsv(out_tx, out_file_tx)

