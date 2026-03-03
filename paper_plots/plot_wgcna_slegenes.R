library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(tidytext)
 
sle_genes <- 
    read_tsv("./data/GCST003155-gwas-credible-sets-studies.tsv") |>
    mutate(topL2G = recode(topL2G, "BLTP3A" = "UHRF1BP1"))
   
gene_annot <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.primary_assembly.annotation.gtf.gz" |>
    rtracklayer::import(feature.type = "gene") |>
    as.data.frame() |>
    as_tibble() |>
    select(gene_id, gene_name)

module_sizes <- 
    read_tsv("../bcell_lowinput/wgcna/data/DN2_modules.tsv") |>
    count(module) |>
    filter(module != "grey") |>
    arrange(desc(n)) |>
    mutate(module_ix = glue::glue("M{1:n()}")) |>
    select(module, module_ix)
    
kme_all_df <- 
    read_tsv("../bcell_lowinput/wgcna/data/DN2_kme.tsv") |>
    select(-grey) |>
    pivot_longer(-gene_id, names_to = "module", values_to = "kme") |>
    left_join(module_sizes) |>
    select(gene_id, module = module_ix, kme)

sle_genes_cormatrix <- 
    gene_annot |>
    filter(gene_name %in% sle_genes$topL2G) |>
    left_join(kme_all_df, join_by(gene_id)) |>
    filter(!is.na(module)) |>
    select(-gene_id) |>
    pivot_wider(names_from = module, values_from = kme) |>
    column_to_rownames("gene_name") |>
    select(all_of(module_sizes$module_ix)) |>
    data.matrix()

png("./sfig_wgcna_sle_genes.png", units = "in", height = 6, width = 4, res = 300)
pheatmap(sle_genes_cormatrix,
	 fontsize = 9, angle_col = 0, cluster_rows = TRUE, cluster_cols = FALSE, 
	 color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(100),
	 display_numbers = round(sle_genes_cormatrix, 2),
	 legend = FALSE)
dev.off()

sle_genes_cormatrix |>
    as_tibble(rownames = "gene") |>
    pivot_longer(-gene, names_to = "module", values_to = "kme") |>
    mutate(kme = round(kme, 2)) |>
    filter(kme > .85) |>
    count(module) |>
    arrange(desc(n))
