library(tidyverse)

##magma <- read.table("./results/Bentham.genes.out", header = TRUE) %>%
##    as_tibble()
#
#magma <- read.table("./results/SleRiskGenes.genes.out", header = TRUE) %>%
#    as_tibble()
#
#sle <- magma %>%
#    select(GENE, SLE = ZSTAT) %>%
#    arrange(desc(SLE))
#
##write_tsv(sle, "./results/Bentham_zscore.tsv")
#write_tsv(sle, "./results/SleRiskGenes_zscore.tsv")
#
# MAGMA scores from scDRS paper
magma <- read.table("./data/MAGMA_v108_GENE_10_ZSTAT_for_scDRS.txt") %>%
    rownames_to_column("gene_name") %>%
    as_tibble() %>%
    select(gene_name, SLE = "PASS_Lupus") %>%
    drop_na() %>%
    arrange(desc(SLE))

gencode <- read_tsv("../../data/gencode.v19.genes.tsv") %>%
    mutate(gene_id = sub("\\.\\d+$", "", gene_id)) %>%
    select(gene_id, gene_name)

gencode38 <- read_tsv("./data/gencode_genes/gencode.v38.GRCh38.bed", col_names = FALSE) %>%
    extract(X4, c("gene_id", "gene_name"), "(ENSG[0-9.]+)_(.+)") %>%
    mutate(gene_id = sub("\\.\\d+$", "", gene_id)) %>%
    select(gene_id, gene_name)
    
gene_ids <- 
    bind_rows(gencode, gencode38) %>%
    distinct() %>%
    add_count(gene_name) %>%
    filter(n == 1) %>%
    select(-n)

out <- left_join(magma, gene_ids) %>% 
    drop_na() %>%
    select(GENE = gene_id, SLE)
    
write_tsv(out, "./data/Bentham_zscore_from_scDRS.tsv")
write_tsv(gene_ids, "./data/gencode_genes/gene_ids.tsv")
