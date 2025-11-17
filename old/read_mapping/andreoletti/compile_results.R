library(tidyverse)

annotations <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.txToGene.tsv" %>%
    read_tsv()

sample_ids <- 
    "/lab-share/IM-Gutierrez-e2/Public/External_datasets/Andreoletti/SRR_Acc_List_Bcells.txt" %>%
    read_lines()

quant_df <- file.path("results", sample_ids, "quant.sf") %>%
    setNames(sample_ids) %>%
    map_dfr(read_tsv, .id = "sampleid")

quant_genes_df <- quant_df %>%
    left_join(annotations, by = c("Name" = "transcript_id")) %>%
    group_by(sampleid, gene_id, gene_name) %>%
    summarise(tpm = sum(TPM)) %>%
    ungroup()

write_tsv(quant_genes_df, "compiled_expression.tsv")
