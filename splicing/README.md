README
================

## Packages

``` r
library(tidyverse)
library(tidytext)
library(ggsci)
library(ggrepel)
library(cowplot)
```

## Functions

``` r
read_leaf <- function(cell_type) {
    
    sig <- 
        file.path("./results", cell_type, "leafcutter_cluster_significance.txt") %>%
        read_tsv() %>%
        separate(cluster, c("chr", "cluster"), sep = ":") %>%
        select(cluster, p.adjust, genes)

    eff <- 
        file.path("./results", cell_type, "leafcutter_effect_sizes.txt") %>%
        read_tsv() %>%
        mutate(cluster = sub("^.*(clu.*)$", "\\1", intron)) %>%
        select(cluster, logef, deltapsi)

    inner_join(sig, eff)
}
```

``` r
cell_types <- c("rN", "aN", "T3", "SM", "DN") %>%
    setNames(., .)

leaf_df <- map_df(cell_types, read_leaf, .id = "cell_type") %>%
    group_by(cell_type, cluster) %>%
    slice(which.max(abs(deltapsi))) %>%
    ungroup()

dtu <- 
    file.path("./results", cell_types, "dtu_perTx.tsv") %>%
    setNames(cell_types) %>%
    map_df(~read_tsv(.) %>%
               select(gene_name, tx_id = feature_id, adj_pvalue),
           .id = "cell_type")
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## GO enrichment

GO enrichment analysis on genes with at least one cluster with absolute
delta PSI &gt; 0.1.

``` r
go <- read_tsv("./results/enrichment.tsv") %>%
    arrange(Description)
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
