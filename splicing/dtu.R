library(tidyverse)
library(tximport)
library(DRIMSeq)

# Cell type
cell_type <- commandArgs(TRUE)[1]
#cell_type <- "DN"

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
		   min_samps_feature_expr = 5, min_feature_expr = 10,
		   min_samps_feature_prop = 5, min_feature_prop = 0.1,
		   min_samps_gene_expr = 10, min_gene_expr = 10)

table(table(counts(d_filt)$gene_id))

## DEXSeq
library(DEXSeq)

sample_data <- DRIMSeq::samples(d_filt)

count_data <- round(as.matrix(counts(d_filt)[, -c(1,2)]))

dxd <- DEXSeqDataSet(countData = count_data,
		     sampleData = sample_data,
		     design = ~sample + exon + condition:exon,
		     featureID = counts(d_filt)$feature_id,
		     groupID = counts(d_filt)$gene_id)

colData(dxd)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet = TRUE)
dxd <- testForDEU(dxd, reducedModel = ~sample + exon)

dxr <- DEXSeqResults(dxd, independentFiltering = FALSE)
qval <- perGeneQValue(dxr)
dxr_g <- data.frame(gene = names(qval), qval)

## stageR
library(stageR)

strp <- function(x) substr(x,1,15)
pConfirmation <- matrix(dxr$pvalue,ncol=1)
dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")
pScreen <- qval
names(pScreen) <- strp(names(pScreen))
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")]) %>%
    mutate_at(vars(1:2), strp) 

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                 onlySignificantGenes=TRUE)
})

write_tsv(dex.padj, file.path("./results", cell_type, "dtu.tsv"))
