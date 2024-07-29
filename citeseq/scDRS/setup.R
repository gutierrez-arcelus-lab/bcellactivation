library(tidyverse)

# GWAS ########################################################################
if (!file.exists("data/magma_input")) dir.create("data/magma_input")
if (!file.exists("output")) dir.create("output")

bentham_stats <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/harmonized/26502338-GCST003156-EFO_0002690-build37.f.tsv.gz" |>
    data.table::fread() |>
    as_tibble() |>
    select(chrom = chromosome, pos = base_pair_location, rsid = variant_id, p = p_value) |>
    filter(chrom != "Y",
	   !(chrom == 6 & between(pos, 25e6, 34e6))) |>
    mutate(chrom = recode(chrom, "X" = "23"),
	   chrom = as.integer(chrom))

bentham_stats |>
    select(rsid, chrom, pos) |>
    write_tsv("./data/magma_input/sle_bentham_snp_loc.tsv", col_names = FALSE)

bentham_stats |>
    select(SNP = rsid, P = p) |>
    write_tsv("./data/magma_input/sle_bentham_snp_pval.tsv")

# Gene locations from Gencode ################################################
gencode <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v19.chr_patch_hapl_scaff.annotation.gtf" |>
    read_tsv(comment = "##", col_names = FALSE) |>
    filter(X3 == "gene", X1 %in% paste0("chr", c(1:22, "X"))) |>
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_id = str_remove(gene_id, "\\.\\d+$")) |>
    select(gene_id, chrom = X1, start = X4, end = X5, strand = X7) |>
    mutate(chrom = str_remove(chrom, "^chr"),
	   chrom = recode(chrom, "X" = "23"))

write_tsv(gencode, "./data/magma_input/gene_loc_hg19.tsv", col_names = FALSE)


# Create H5 object from Seurat ################################################
if (!file.exists("data/scdrs_input")) dir.create("data/scdrs_input")

library(Seurat)
library(SeuratData)
library(SeuratDisk)

seurat <- readRDS("../data/seurat_qced.rds")
seurat <- DietSeurat(seurat, dimreducs = NULL, assays = "RNA", counts = FALSE)

seurat@commands <- list()
seurat@assays$RNA@var.features <- list()

SaveH5Seurat(seurat, filename = "./data/scdrs_input/bcells.h5Seurat")
Convert("./data/scdrs_input/bcells.h5Seurat", dest = "h5ad")
