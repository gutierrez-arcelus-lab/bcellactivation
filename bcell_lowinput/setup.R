library(tidyverse)

meta <- 
    "/lab-share/IM-Gutierrez-e2/Public/Lab_datasets/B_cells_lowinput/metadata/description_hsapiens_qc.tsv" |>
    read_tsv() |>
    mutate(sample_id = str_remove(sample_id, "\\.rep\\.\\d$")) |>
    group_by(sample_id, stim, time) |>
    summarise_at(vars(fq1, fq2), ~paste(., collapse = ",")) |>
    ungroup() |>
    filter(sample_id != "BLANK") |>
    unite("id", c(sample_id, stim, time), sep = "_")

write_tsv(meta, "./data/metadata_pooledreps.tsv", col_names = FALSE)

# Curated list of SLE genes
selected_genes <-
    c(
      "ARID5B",
      "ATG5",
      "BANK1",
      "BLK", 
      "CD44",
      "CSK",
      "DHCR7",
      "ETS1",
      "FCGR2A",
      "LYST",
      "IFIH1",
      "IKBKE",
      "IKZF1",
      "IKZF2",
      "IKZF3",
      "IL10",
      "IL12A",
      "IRAK1",
      "IRF5",
      "IRF7",
      "IRF8",
      "ITGAM",
      "ITGAX",
      "JAZF1",
      "MECP2",
      "MIR146A",
      "MIR3142HG",
      "NCF2",
      "PTPN22",
      "PRDM1",
      "PXK",
      "RAD51B",
      "SH2B3",
      "SLC15A4",
      "SOCS1",
      "SPRED2",
      "STAT1", 
      "STAT4",
      "TASL",
      "TCF7",
      "TNFAIP3",
      "TNFSF4",
      "TNIP1",
      "TREX1",
      "TYK2",
      "UHRF1BP1",
      "UBE2L3",
      "WDFY4"
    )


write_lines(selected_genes, "./data/sle_curated_genes.txt")
