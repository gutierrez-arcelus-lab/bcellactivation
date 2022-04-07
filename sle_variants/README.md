Sle GWAS variants
================

### Packages

``` r
library(tidyverse)
library(readxl)
library(rvest)
```

### Functions

``` r
parse_pvalues <- function(p) {

    p %>%
        gsub("\\s", "", .) %>%
        sub("^([0-9.]+)[^0-9.]10[^0-9.](\\d+)$", "\\1e-\\2", .) %>%
        parse_number()
}

read_langefeld <- function(tab_path, sheet_number) {
    read_excel(tab_path, sheet_number, skip = 2) %>%
    select(snp_id = 1, 
           chr = 2, 
           pos = 3, 
           gene_region = 4, 
           region_rank = 5, 
           ref_allele = 6, 
           p = `P-value`,
           or = `OR (95% CI)`,
           p_stepwise = `Regional Stepwise  P-value`) %>% 
    filter(!is.na(pos)) %>%
    mutate(snp_id = gsub("\\s|[a-zA-Z,]+$", "", snp_id),  
           p = parse_pvalues(p),
           p_stepwise = parse_pvalues(p_stepwise),
           or = parse_number(or)) %>%
    select(snp_id, chr, pos, gene_region, region_rank, ref_allele, p, p_stepwise, or)
}

read_bentham <- function(url_address) {

    read_html(url_address) %>%
    html_node("table") %>%
    html_table(header = TRUE, fill = TRUE) %>%
    select(snp_id = 1, chr = 2, pos = 3, locus = 4, p = 11, or = 12) %>%
    slice(-1) %>%
    as_tibble() %>%
    mutate(snp_id = gsub("\\s|[a-zA-Z,]+$", "", snp_id),
           pos = parse_number(pos),
           p = parse_pvalues(p),
           or = parse_number(or))
}
```

### Import Langefeld et al. data

``` r
tier1_ea <- read_langefeld("./paper_data/langefeld_tableS2.xlxs", 1)
tier2_ea <- read_langefeld("./paper_data/langefeld_tableS2.xlxs", 2)
tier3_ea <- read_langefeld("./paper_data/langefeld_tableS2.xlxs", 3)

langefeld_df <- bind_rows("1" = tier1_ea, "2" = tier2_ea, "3" = tier3_ea, .id = "tier") %>%
    arrange(tier, region_rank, p) %>%
    mutate(chr = sub("^(\\d+)[qp].*$", "\\1", chr),
           chr = factor(chr, levels = c(1:22, "X")))

langefeld_df
```

    # A tibble: 798 x 10
       tier  snp_id     chr         pos gene_region region_rank ref_allele        p p_stepwise    or
       <chr> <chr>      <fct>     <dbl> <chr>             <dbl> <chr>         <dbl>      <dbl> <dbl>
     1 1     rs12706861 7     128616582 IRF5-TNPO3            1 T          3.85e-71  NA         1.76
     2 1     rs35000415 7     128585616 IRF5-TNPO3            1 T          5.67e-70   2.98e-35  1.73
     3 1     rs3807307  7     128579202 IRF5-TNPO3            1 C          3.75e-62   1.33e-23  1.46
     4 1     rs7808907  7     128584084 IRF5-TNPO3            1 T          1.44e-31   8.33e- 6  0.77
     5 1     rs7582694  2     191970120 STAT4                 2 C          4.3 e-69  NA         1.56
     6 1     rs7568275  2     191966452 STAT4                 2 C          4.51e-68   7.67e-70  1.55
     7 1     rs6715106  2     191913034 STAT4                 2 G          8.33e-15   1.17e- 7  0.67
     8 1     rs932169   2     191929278 STAT4                 2 G          3.29e- 6   2.12e-14  1.2 
     9 1     rs1143679  16     31276811 ITGAM                 3 A          2.78e-62  NA         1.72
    10 1     rs34572943 16     31272353 ITGAM                 3 A          2.63e-58   2.32e-62  1.67
    # … with 788 more rows

### Import Bentham et al. data

``` r
bentham_df <- read_bentham("https://www.nature.com/articles/ng.3434/tables/1") %>%
    mutate(chr = factor(chr, levels = c(1:22, "X")))

bentham_df
```

    # A tibble: 43 x 6
       snp_id     chr         pos locus              p    or
       <chr>      <fct>     <dbl> <chr>          <dbl> <dbl>
     1 rs2476601  1     114377568 PTPN22     1.1 e- 28  1.43
     2 rs1801274  1     161479745 FCGR2A     1.04e- 12  1.16
     3 rs704840   1     173226195 TNFSF4     3.12e- 19  1.22
     4 rs17849501 1     183542323 SMG7, NCF2 3.45e- 88  2.1 
     5 rs3024505  1     206939904 IL10       4.64e-  9  1.17
     6 rs9782955  1     236039877 LYST       1.25e-  9  1.16
     7 rs6740462  2      65667272 SPRED2     2.67e-  5  1.1 
     8 rs2111485  2     163110536 IFIH1      1.27e- 11  1.15
     9 rs11889341 2     191943742 STAT4      5.59e-122  1.73
    10 rs3768792  2     213871709 IKZF2      1.21e- 13  1.24
    # … with 33 more rows

### Getting GRCh38 positions and p-values from the summary statistics

``` r
bentham_stats <- read_tsv("./summ_stats/bentham_GRCh38.tsv.gz")
```

``` r
bentham_38 <- bentham_stats %>%
    select(chr = chromosome, pos = base_pair_location, snp_id = variant_id, p_value) %>%
    mutate(chr = ifelse(chr == 23, "X", chr),
           chr = factor(chr, levels = c(1:22, "X"))) %>%
    left_join(bentham_df, ., by = c("chr", "snp_id")) %>%
    select(snp_id, chr, locus, pos = pos.y, p_value)
    
bentham_38
```

    # A tibble: 43 x 5
       snp_id     chr   locus            pos  p_value
       <chr>      <fct> <chr>          <dbl>    <dbl>
     1 rs2476601  1     PTPN22     113834946 8.38e-13
     2 rs1801274  1     FCGR2A     161509955 5.56e-11
     3 rs704840   1     TNFSF4     173257056 1.41e-13
     4 rs17849501 1     SMG7, NCF2 183573188 1.81e-59
     5 rs3024505  1     IL10       206766559 2.21e- 3
     6 rs9782955  1     LYST       235876577 7.92e- 4
     7 rs6740462  2     SPRED2      65440138 1.76e- 8
     8 rs2111485  2     IFIH1      162254026 3.69e- 6
     9 rs11889341 2     STAT4      191079016 1.12e-65
    10 rs3768792  2     IKZF2      213006985 3.78e- 8
    # … with 33 more rows
