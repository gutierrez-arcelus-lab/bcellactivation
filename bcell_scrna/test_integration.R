library(tidyverse)
library(cowplot)
library(scico)
library(Seurat)


pilot1 <- readRDS("./data/pilot1.Rdata")
pilot2 <- readRDS("./data/pilot2.Rdata")

pilot1 <- pilot1 %>%
    FindVariableFeatures(nfeatures = 2000, selection.method = "vst") %>%
    ScaleData(., features = rownames(.)) %>%
    RunPCA(., features = VariableFeatures(.)) %>%
    RunUMAP(dims = 1:20, verbose = FALSE)

pilot2 <- pilot2 %>%
    FindVariableFeatures(nfeatures = 2000, selection.method = "vst") %>%
    ScaleData(., features = rownames(.)) %>%
    RunPCA(., features = VariableFeatures(.)) %>%
    RunUMAP(dims = 1:20, verbose = FALSE)

umap_1 <- as.data.frame(pilot1@reductions$umap@cell.embeddings) %>%
    as_tibble(rownames = "barcode")

umap_2 <- as.data.frame(pilot2@reductions$umap@cell.embeddings) %>%
    as_tibble(rownames = "barcode")

adt_1 <- pilot1@assays$ADT@data %>%
    .[grep("ctrl", rownames(.), ignore.case = TRUE, value = TRUE), ] %>%
    as_tibble(rownames = "ab") %>%
    pivot_longer(-ab, names_to = "barcode") %>%
    group_by(ab) %>%
    mutate(q01 = quantile(value, 0.01),
	   q99 = quantile(value, 0.99),
	   value = case_when(value < q01 ~ q01,
			     value > q99 ~ q99,
			     TRUE ~ value)) %>%
    ungroup() %>%
    select(ab, barcode, value)

adt_2 <- pilot2@assays$ADT@data %>%
    .[grep("ctrl", rownames(.), ignore.case = TRUE, value = TRUE), ] %>%
    as_tibble(rownames = "ab") %>%
    pivot_longer(-ab, names_to = "barcode") %>%
    group_by(ab) %>%
    mutate(q01 = quantile(value, 0.01),
	   q99 = quantile(value, 0.99),
	   value = case_when(value < q01 ~ q01,
			     value > q99 ~ q99,
			     TRUE ~ value)) %>%
    ungroup() %>%
    select(ab, barcode, value)

plot_1 <- umap_1 %>%
    left_join(adt_1, by = "barcode") %>%
    ggplot(aes(UMAP_1, UMAP_2)) +
	geom_point(aes(color = value), size = .25) +
	scale_color_scico("", palette = "lajolla",
			  guide = guide_colorbar(barwidth = .5),
			  limits = c(0, max(c(adt_1$value, adt_2$value))),
			  breaks = scales::pretty_breaks(6)) +
	theme_bw() +
	theme(panel.grid = element_blank()) +
	facet_wrap(~ab) +
	labs(title = "Pilot 1")

plot_2 <- umap_2 %>%
    left_join(adt_2, by = "barcode") %>%
    ggplot(aes(UMAP_1, UMAP_2)) +
	geom_point(aes(color = value), size = .25) +
	scale_color_scico("", palette = "lajolla",
			  guide = guide_colorbar(barwidth = .5),
			  limits = c(0, max(c(adt_1$value, adt_2$value))),
			  breaks = scales::pretty_breaks(6)) +
	theme_bw() +
	theme(panel.grid = element_blank()) +
	facet_wrap(~ab) +
	labs(title = "Pilot 2")

out <- plot_grid(plot_1, plot_2, ncol = 1)


ggsave("./mouseadts.png", out, width = 8, height = 10) 

# Integrate data

pilot_list <- list("pilot 1" = pilot1, "pilot 2" = pilot2)
  
features <- SelectIntegrationFeatures(object.list = pilot_list)

pilot_anchors <- FindIntegrationAnchors(object.list = pilot_list,
					anchor.features = features,
                                        reduction = "cca",
                                        k.anchor = 20,
                                        dims = 1:20)

pilots_integrated <- 
    IntegrateData(anchorset = pilot_anchors, 
		  dims = 1:20,
		  features.to.integrate = rownames(pilot1@assays$RNA@counts))

pilots_integrated <- pilots_integrated %>%
    ScaleData(., features = rownames(.), verbose = FALSE) %>%
    RunPCA(., features = VariableFeatures(.), verbose = FALSE) %>%
    RunUMAP(dims = 1:20, verbose = FALSE, seed.use = 1) %>%
    FindNeighbors(dims = 1:20, verbose = FALSE) %>%
    FindClusters(resolution = 0.5, verbose = FALSE)

meta_df <- pilots_integrated@meta.data %>% 
    as_tibble(rownames = "barcode") %>%
    select(barcode, pilot = orig.ident, stim = HTO_maxID)

umap_int <- pilots_integrated@reductions$umap@cell.embeddings %>%
    as_tibble(rownames = "barcode") %>%
    left_join(meta_df)

pilot1_stim <- meta_df %>%
    filter(pilot == "pilot1") %>%
    pull(stim) %>%
    unique() %>%
    str_sort(numeric = TRUE) %>%
    fct_relevel("Res 0hr", after = 0) %>%
    fct_relevel("Res 24hr", after = 1) %>%
    levels()

pilot2_stim <- meta_df %>%
    filter(pilot == "pilot2") %>%
    pull(stim) %>%
    unique() %>%
    str_sort(numeric = TRUE) %>%
    fct_relevel("day0", after = 0) %>%
    fct_relevel("IL4 24hr", after = 1) %>%
    fct_relevel("TLR7 24hr", after = 4) %>%
    fct_relevel("TLR7 72hr", after = 5) %>%
    levels()

library(RColorBrewer)

pilot1_pal <- c("grey", "slategrey", brewer.pal(n = 4, "Paired")) %>%
    setNames(pilot1_stim)

pilot2_pal <- c("grey", "slategrey", brewer.pal(n = 12, "Paired")) %>%
    setNames(pilot2_stim)

int_p1 <- 
    ggplot(filter(umap_int, pilot == "pilot1"),
	   aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = stim), size = .5) +
    scale_color_manual(values = pilot1_pal) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    labs(title = "Pilot 1") +
    guides(color = guide_legend(override.aes = list(size = 2)))

int_p2 <- 
    ggplot(filter(umap_int, pilot == "pilot2"),
	   aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = stim), size = .5) +
    scale_color_manual(values = pilot2_pal) +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    labs(title = "Pilot 2") +
    guides(color = guide_legend(override.aes = list(size = 2)))

out_int <- plot_grid(int_p1, int_p2, nrow = 1) + 
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./integrated_pilots.png", out_int, width = 10, height = 4) 



