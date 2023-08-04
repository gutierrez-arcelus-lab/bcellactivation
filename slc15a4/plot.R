library(tidyverse)
library(ggsci)

coloc_df <- read_tsv("./data/coloc_results.tsv")

plot_df <- 
    coloc_df |>
    group_by(study, tissue, gene_id, mol_phenotype) |>
    slice_max(h4) |>
    ungroup() |>
    mutate(study_id = sprintf("%s (%s)", study, tissue),
	   study_id = gsub("_", " ", study_id)) |>
    select(study_id, mol_phenotype, gene_id, gene_name, pp4 = h4) |>
    group_by(gene_id) |>
    filter(any(pp4 > .2)) |>
    group_by(study_id) |>
    filter(any(pp4 > .2)) |>
    ungroup() |>
    mutate_at(vars(mol_phenotype), ~factor(., levels = c("ge", "exon", "microarray", "tx", "txrev", "leafcutter"))) |>
    arrange(desc(pp4)) |>
    mutate(study_id = fct_inorder(study_id),
	   gene_name = fct_relevel(gene_name, "SLC15A4", after = 0))


p <- ggplot(plot_df, aes(x = pp4, y = study_id)) +
    geom_vline(xintercept = .9, linetype = 2, linewidth = .2, alpha = .5) +
    geom_point(aes(color = gene_name), size = 3) +
    facet_wrap(~mol_phenotype, nrow = 1) +
    scale_x_continuous(limits = c(0, 1),
		       breaks = seq(0, 1),
		       labels = c("0", "1")) +
    scale_color_npg() +
    theme(panel.background = element_rect(fill = "grey96")) +
    labs(x = "Posterior probability of colocalization", 
	 y = NULL,
	 color = "Gene")

ggsave("./plots/colocs.png", p, width = 10, height = 7, dpi = 600)

eqtl_cat <- read_tsv("./data/eqtl_catalogue_associations.tsv.gz")

eqtl_cat |> 
    filter(dataset_id == "QTD000089", grepl("clu_21371", molecular_trait_id)) |>
    arrange(pvalue) |>
    print(width = Inf)

eqtl_cat |> 
    filter(dataset_id == "QTD000089", grepl("clu_21371", molecular_trait_id)) |>
    filter(position == 128788426)



