library(tidyverse)

#genes <- read_tsv("./data/ncbi_genes/NCBI37.3.gene.loc", col_names=FALSE)
#genes <- read_tsv("./data/gencode_genes/genes.v38.hg19.loc", col_names = FALSE)

magma <- read.table("./results/Bentham.genes.out", header = TRUE) %>%
    as_tibble()

sle <- magma %>%
    #left_join(select(genes, X1, X6), by = c("GENE" = "X1")) %>%
    select(GENE, SLE = ZSTAT) %>%
    arrange(desc(SLE))

write_tsv(sle, "./results/Bentham_zscore.tsv")
