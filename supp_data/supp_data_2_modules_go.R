# ==============================================================================
# Description:  Aggregates Gene Ontology (GO) enrichment results across all WGCNA 
#               networks into a single supplementary Excel file. 
# Input:        1. *_modules.tsv (Gene module assignments per condition)
#               2. *_go.tsv (GO enrichment results per condition)
# Output:       supplementary_data_2.xlsx (Multi-sheet Excel workbook)
# ==============================================================================

library(tidyverse)
library(glue)
library(writexl)

# Define the standard stimulation conditions to iterate over
my_stims <- c("CD40L", "TLR9", "TLR7", "BCR", "BCR-TLR7", "DN2")

# ------------------------------------------------------------------------------
# 1. Calculate Module Sizes and Generate Module IDs
# ------------------------------------------------------------------------------
module_sizes <-
    glue("../01_rnaseq_lowinput/3_wgcna/data/{my_stims}_modules.tsv") |>
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

# ------------------------------------------------------------------------------
# 2. Import Gene Ontology Results
# ------------------------------------------------------------------------------
go_all <-
    glue("../01_rnaseq_lowinput/3_wgcna/data/{my_stims}_go.tsv") |>
    setNames(my_stims) |>
    map_dfr(read_tsv, .id = "stim") |>
    select(stim, module, Description, geneID, pvalue, p.adjust, qvalue)

# ------------------------------------------------------------------------------
# 3. Merge and Format Data
# ------------------------------------------------------------------------------
stab <- 
    go_all |>
    left_join(module_sizes, join_by(stim, module)) |>
    mutate(stim = factor(stim, levels = my_stims)) |>
    arrange(stim, ix) |>
    select(stim, module, ix, Description, geneID, pvalue, p.adjust, qvalue)

# ------------------------------------------------------------------------------
# 4. Export Multi-Sheet Excel Workbook
# ------------------------------------------------------------------------------
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

# Passing a named list to write_xlsx automatically generates an Excel workbook 
# where each list element becomes its own named sheet.
write_xlsx(stab_list, "./data/supplementary_data_2_modules_go.xlsx")
