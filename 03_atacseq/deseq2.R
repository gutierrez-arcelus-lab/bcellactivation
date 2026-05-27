# Code borrowed from the featurecounts_deseq2.r script in 
# nf-core atacseq pipeline v1.2.2

library(DESeq2)
library(vsn)
library(BiocParallel)
library(tidyverse)
library(glue)

# meta data
samplesheet <- 
    read_csv("./samplesheet.csv") |>
    mutate(f1 = basename(fastq_1),
	   donor_id = str_extract(f1, "20221025_([^_]+).*", group = 1),
	   replic = paste0("REP", replicate)) |>
    select(condition = sample, replic, donor_id) |>
    distinct()

write_tsv(samplesheet, "./results_deseq2/metadata.tsv")


# Count data
dat <- 
    "./results/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.featureCounts.txt" |>
    read_table(skip = 1)

colnames(dat) <- gsub("\\.mLb\\.\\clN\\.sorted\\.bam", "", colnames(dat))

# Interval info
interval_table <-
    dat |>
    select(Geneid:Length) |>
    column_to_rownames("Geneid")

# Counts
samples_vec <- 
    c(
    sort(grep("unst_", colnames(dat), value = TRUE)),
    sort(grep("IL4_", colnames(dat), value = TRUE)),
    sort(grep("TLR7_", colnames(dat), value = TRUE)),
    sort(grep("BCR_", colnames(dat), value = TRUE)),
    sort(grep("DN2_", colnames(dat), value = TRUE))
    )

# Remove pooled sample of 3 donors
pool_samples <- 
    samplesheet |>
    filter(donor_id == "3donors") |>
    unite("sample_id", c(condition, replic), sep = "_") |>
    pull(sample_id)

count_table <- 
    dat |>
    select(Geneid, all_of(samples_vec)) |>
    select(-all_of(pool_samples)) |>
    column_to_rownames("Geneid")

## RUN DESEQ2
coldata <- 
    tibble(sample_id = colnames(count_table)) |>
    separate(sample_id, c("stim", "tp", "replic"), sep = "_", remove = FALSE) |>
    unite("condition", c(stim, tp), sep = "_") |>
    left_join(samplesheet) |>
    select(-replic) |>
    mutate_at(vars(condition, donor_id), fct_inorder) |>
    column_to_rownames("sample_id")

dds <- 
    DESeqDataSetFromMatrix(countData = round(count_table), 
			   colData = coldata, 
			   design = ~ condition + donor_id)

dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(4))
rld <- rlog(dds)

save(dds, rld, file = "./results_deseq2/deseq2_data.Rdata")


## PCA
walk(c(500, 1000, 2500, 5000, 10000), 
     function(x) {
	 DESeq2::plotPCA(rld, intgroup = "condition", ntop = x, returnData = TRUE) |>
	 write_rds(glue("./results_deseq2/pcadata_{x}peaks.rds"))
     })

# Differential accessibility analysis
comparisons <- 
    coldata$condition |>
    unique() |>
    as.character() |>
    combn(2)

for (idx in 1:ncol(comparisons)) {

    control_group <- comparisons[1, idx]
    treat_group <- comparisons[2, idx]
    control_samples <- filter(coldata, condition == control_group) |> rownames()
    treat_samples <- filter(coldata, condition == treat_group) |> rownames()

    results_da <- 
	results(dds, contrast = c("condition", treat_group, control_group))

    results_out <- 
	cbind(interval_table, as.data.frame(results_da)) |>
	as_tibble(rownames = "interval")

    write_tsv(results_out, glue("./results_deseq2/{treat_group}vs{control_group}.tsv"))

}
