library(tidyverse)
library(WGCNA)

counts_dn2 <- read_rds("./data/DN2_counts.rds")
counts_bcr <- read_rds("./data/BCR_counts.rds")

modules_dn2 <-
    read_tsv("./data/DN2_modules.tsv") |>
    deframe() |>
    {function(x) x[colnames(counts_dn2)]}()

modules_bcr <-
    read_tsv("./data/BCR_modules.tsv") |>
    deframe() |>
    {function(x) x[colnames(counts_bcr)]}()

multi_expr <- list()
multi_expr[[1]] <- list(data = counts_dn2)
multi_expr[[2]] <- list(data = counts_bcr)
names(multi_expr) <- c("DN2", "BCR")

color_list <- list("DN2" = modules_dn2, "BCR" = modules_bcr)

mp <- 
    modulePreservation(multi_expr, color_list, 
		       referenceNetworks = 1,
		       testNetworks = 2,
		       loadPermutedStatistics = FALSE,
		       nPermutations = 100,
		       verbose = 3)

str(mp)

# Correlation of module eigen genes
library(RColorBrewer)
library(pheatmap)

eigen_dn2 <- 
    read_tsv("./data/DN2_eigen.tsv") |>
    select(-MEgrey) |>
    separate(sample_name, c("donor_id", "stim", "time"), sep = "_") |>
    unite("sample_id", c(donor_id, time), sep = "_") |>
    select(-stim) |>
    pivot_longer(-sample_id, names_to = "module_dn2", values_to = "values_dn2")
    
eigen_bcr <- 
    read_tsv("./data/BCR_eigen.tsv") |>
    select(-MEgrey) |>
    separate(sample_name, c("donor_id", "stim", "time"), sep = "_") |>
    unite("sample_id", c(donor_id, time), sep = "_") |>
    select(-stim) |>
    pivot_longer(-sample_id, names_to = "module_bcr", values_to = "values_bcr")

eigen_tlr7 <- 
    read_tsv("./data/TLR7_eigen.tsv") |>
    select(-MEgrey) |>
    separate(sample_name, c("donor_id", "stim", "time"), sep = "_") |>
    unite("sample_id", c(donor_id, time), sep = "_") |>
    select(-stim) |>
    pivot_longer(-sample_id, names_to = "module_tlr7", values_to = "values_tlr7")


# DN2 vs BCR
mod_cors <-
    inner_join(eigen_dn2, eigen_bcr, 
	       join_by(sample_id),
	       relationship = "many-to-many") |>
    group_by(module_dn2, module_bcr) |>
    summarise(r = cor(x = values_dn2, y = values_bcr)[,1]) |>
    ungroup() |>
    mutate_at(vars(module_dn2, module_bcr), ~str_remove(., "ME")) |>
    pivot_wider(names_from = module_bcr, values_from = r) |>
    column_to_rownames("module_dn2")


png("./plots/heatmap_modules_DN2_BCR.png", 
    units = "in", height = 5, width = 5, res = 300)
pheatmap(mod_cors,
	 fontsize = 9, angle_col = 90,
	 display_numbers = TRUE,
	 color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),
	 legend_breaks = seq(-1, 1, by = .2))
dev.off()

mod_cors_long <- 
    mod_cors |>
    rownames_to_column("module_dn2") |>
    pivot_longer(-module_dn2, names_to = "module_bcr") |>
    group_by(module_dn2) |>
    slice_max(value) |>
    ungroup()

bcr_kme <- 
    read_tsv("./data/BCR_kme.tsv") |>
    pivot_longer(-gene_id, names_to = "module", values_to = "kme") |>
    filter(kme >= 0.9) |>
    mutate(module = str_remove(module, "kME")) |>
    group_by(gene_id) |>
    slice_max(kme) |>
    ungroup() |>
    select(gene_id, module)

dn2_kme <- 
    read_tsv("./data/DN2_kme.tsv") |>
    pivot_longer(-gene_id, names_to = "module", values_to = "kme") |>
    filter(kme >= 0.9) |>
    mutate(module = str_remove(module, "kME")) |>
    group_by(gene_id) |>
    slice_max(kme) |>
    ungroup() |>
    select(gene_id, module)

inner_join(dn2_kme, bcr_kme, join_by(gene_id), suffix = c("_dn2", "_bcr")) |>
    count(module_dn2, module_bcr, sort = TRUE) |>
    print(n = Inf)

gene_ids <- 
    "../results/edger/results.tsv" |>
    read_tsv() |>
    distinct(gene_id, gene_name)

inner_join(dn2_kme, bcr_kme, join_by(gene_id), suffix = c("_dn2", "_bcr")) |>
    filter(module_dn2 == "brown", module_bcr == "brown") |>
    left_join(gene_ids)

    
###############################################################################
# BCR vs TLR7

mod_cors_tlr7 <-
    inner_join(eigen_bcr, eigen_tlr7, 
	       join_by(sample_id),
	       relationship = "many-to-many") |>
    group_by(module_bcr, module_tlr7) |>
    summarise(r = cor(x = values_bcr, y = values_tlr7)[,1]) |>
    ungroup() |>
    mutate_at(vars(module_bcr, module_tlr7), ~str_remove(., "ME")) |>
    pivot_wider(names_from = module_tlr7, values_from = r) |>
    column_to_rownames("module_bcr")


png("./plots/heatmap_modules_BCR_TLR7.png", 
    units = "in", height = 5, width = 5, res = 300)
pheatmap(mod_cors_tlr7,
	 fontsize = 9, angle_col = 90,
	 display_numbers = TRUE,
	 color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),
	 legend_breaks = seq(-1, 1, by = .2))
dev.off()

mod_cors_long <- 
    mod_cors |>
    rownames_to_column("module_dn2") |>
    pivot_longer(-module_dn2, names_to = "module_bcr") |>
    group_by(module_dn2) |>
    slice_max(value) |>
    ungroup()

bcr_kme <- 
    read_tsv("./data/BCR_kme.tsv") |>
    pivot_longer(-gene_id, names_to = "module", values_to = "kme") |>
    filter(kme >= 0.9) |>
    mutate(module = str_remove(module, "kME")) |>
    group_by(gene_id) |>
    slice_max(kme) |>
    ungroup() |>
    select(gene_id, module)

dn2_kme <- 
    read_tsv("./data/DN2_kme.tsv") |>
    pivot_longer(-gene_id, names_to = "module", values_to = "kme") |>
    filter(kme >= 0.9) |>
    mutate(module = str_remove(module, "kME")) |>
    group_by(gene_id) |>
    slice_max(kme) |>
    ungroup() |>
    select(gene_id, module)

inner_join(dn2_kme, bcr_kme, join_by(gene_id), suffix = c("_dn2", "_bcr")) |>
    count(module_dn2, module_bcr, sort = TRUE) |>
    print(n = Inf)

gene_ids <- 
    "../results/edger/results.tsv" |>
    read_tsv() |>
    distinct(gene_id, gene_name)

inner_join(dn2_kme, bcr_kme, join_by(gene_id), suffix = c("_dn2", "_bcr")) |>
    filter(module_dn2 == "brown", module_bcr == "brown") |>
    left_join(gene_ids)

