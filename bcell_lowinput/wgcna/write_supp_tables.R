library(tidyverse)
library(glue)
library(writexl)

my_stims <- c("CD40L", "TLR9", "TLR7", "BCR", "BCR-TLR7", "DN2")

module_sizes <-
    glue("./data/{my_stims}_modules.tsv") |>
    setNames(my_stims) |>
    map_dfr(~read_tsv(.) |>
	    count(module) |>
	    arrange(desc(n)) |>
	    filter(module != "grey") |>
	    rowid_to_column("ix") |>
	    select(module, ix) |>
	    mutate(ix = paste("Module", ix)),
	    .id = "stim") |>
    mutate(ix = fct_inorder(ix))

go_all <-
    glue("./data/{my_stims}_go.tsv") |>
    setNames(my_stims) |>
    map_dfr(read_tsv, .id = "stim") |>
    select(stim, module, Description, geneID, pvalue, p.adjust, qvalue)

stab <- 
    go_all |>
    left_join(module_sizes, join_by(stim, module)) |>
    mutate(stim = factor(stim, levels = my_stims)) |>
    arrange(stim, ix) |>
    select(stim, module, ix, Description, geneID, pvalue, p.adjust, qvalue)

stab_list <- 
    stab |>
    mutate(stim = recode(stim, 
			 "CD40L" = "CD40c",
			 "TLR9" = "TLR9c",
			 "TLR7" = "TLR7c",
			 "BCR" = "BCRc",
			 "BCR-TLR7" = "BCR-TLR7c",
			 "DN2" = "DN2c")) |>
    {function(x) split(x, x$stim)}() |>
    map(~select(., module = ix, Description:qvalue))

write_xlsx(stab_list, "./data/supplementary_table_go.xlsx")
