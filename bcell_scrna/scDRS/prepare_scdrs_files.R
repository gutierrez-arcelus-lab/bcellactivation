library(tidyverse)

magma <- read.table("./results/Bentham.genes.out", header = TRUE) %>%
    as_tibble()

sle <- magma %>%
    select(GENE, SLE = ZSTAT) %>%
    arrange(desc(SLE))

write_tsv(sle, "./results/Bentham_zscore.tsv")
