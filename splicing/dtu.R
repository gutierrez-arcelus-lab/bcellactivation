library(tidyverse)
library(tximport)
library(DRIMSeq)

# Sample specs
samps <- "./results/DN/groups_file.txt" %>%
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

