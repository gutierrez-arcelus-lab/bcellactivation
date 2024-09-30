library(tidyverse)

meta_data <- 
    read_tsv("./data/coloc_input/eqtl_catalogue_paths.tsv") |>
    unite("study_id", c(study, quant_method, qtl_group), sep = "_", remove = FALSE) |>
    select(study_id, study, tissue_label, condition_label, quant_method)

res <- 
    sprintf("./results/bentham_region%s.tsv", 1:33) |>
    setNames(1:33) |>
    map_df(~read_tsv(.) |>
	   select(study_id = study, gene_id, gene_name, molecular_trait_id, pp4 = h4) |>
	   mutate(study_id = case_when(grepl("BLUEPRINT", study_id) ~ str_remove(study_id, "(SE|PE)_"),
				       TRUE ~ study_id)) |>
	   left_join(meta_data, join_by(study_id)) |>
	   select(study, tissue = tissue_label, condition = condition_label,
		  gene_id, gene_name, molecular_trait_id, molecular_pheno = quant_method, pp4) |>
	   arrange(study, tissue, condition, molecular_pheno, gene_name, desc(pp4)),
	   .id = "region")

write_tsv(res, "./bentham_colocs.tsv")


