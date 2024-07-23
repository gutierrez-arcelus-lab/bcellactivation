library(tidyverse)

annot <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.primary_assembly.annotation.gtf" |>
    read_tsv(comment = "#", col_names = FALSE)

annot |>
    filter(X3 == "gene") |>
    select(X1, X4, X5) |>
    filter(X1 %in% paste0("chr", c(1:22, "X"))) |>
    distinct() |>
    arrange(X1, X4, X5) |>
    write_tsv("../data/gene_regions.bed", col_names = FALSE)
