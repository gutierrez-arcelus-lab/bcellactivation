library(tidyverse)
library(furrr)

annotations <- 
    "/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/gencode.v38.primary_assembly.annotation.gtf.gz" |>
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c") |>
    mutate(X9 = trimws(X9))

transcript_annots <- annotations %>%
    filter(X3 == "transcript") %>%
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	   transcript_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+")) %>%
    distinct(transcript_id, gene_id, gene_name)

gene_annots <- annotations %>%
    filter(X3 == "gene") %>%
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) %>%
    select(chr = X1, start = X4, end = X5, strand = X7, gene_id, gene_name)

ebv_annots <- 
    "../data/ebv.custom.gtf" %>%
    read_tsv(comment = "track", col_names = FALSE, col_types = "c-cii-c-c") %>%
    mutate(X9 = trimws(X9)) %>%
    filter(X3 == "transcript") %>%
    mutate(gene_id = str_extract(X9, "(?<=gene_id \")[^\"]+"),
	   chr = "ebv") %>%
    select(chr, start = X4, end = X5, strand = X7, gene_id) %>%
    group_by(chr, gene_id, strand) %>%
    summarise(start = min(start), end = max(end)) %>%
    ungroup()

plan(multisession, workers = availableCores())

quant_df <- 
    list.files("./results/salmon", pattern = "quant\\.sf", full.names=TRUE, recursive = TRUE) %>%
    setNames(str_extract(., "ERR\\d+")) %>%
    future_map_dfr(read_tsv, .id = "sampleid")


#### CR2
cr2_ids <- gene_annots %>%
    filter(gene_name == "CR2") %>%
    left_join(transcript_annots) %>%
    pull(transcript_id)

cr2_df <- quant_df %>%
    filter(Name %in% cr2_ids) %>%
    select(sampleid, transcript_id = Name, tpm = TPM)

write_tsv(cr2_df, "./results/salmon/cr2_quants.tsv")
####


human_genes_df <- quant_df %>%
    filter(!grepl("^HHV", Name))

quant_tx_df <- 
    human_genes_df |>
    left_join(transcript_annots, by = c("Name" = "transcript_id"))

quant_tx_df


quant_gene_df <- human_genes_df %>%
    left_join(transcript_annots, by = c("Name" = "transcript_id")) %>%
    group_by(sampleid, gene_id) %>%
    summarise(tpm = sum(TPM)) %>%
    ungroup()

expressed_genes_df <- quant_gene_df %>%
    group_by(gene_id) %>%
    filter(mean(tpm > 0) >= 0.5) %>%
    ungroup()

bed_humans <- expressed_genes_df %>%
    left_join(gene_annots, by = c("gene_id")) %>%
    select(chr, start, end, id = gene_id, gid = gene_name, strd = strand, sampleid, tpm) %>%
    pivot_wider(names_from = sampleid, values_from = tpm) %>%
    filter(chr %in% paste0("chr", c(1:22, "X"))) %>%
    mutate(chr = factor(chr, levels = paste0("chr", c(1:22, "X")))) %>%
    arrange(chr, start)

quant_ebv_df <- quant_df %>%
    filter(grepl("^HHV", Name)) %>%
    separate(Name, c("gene_id", "tmp_id"), sep = ":") %>%
    mutate(gene_id = sub("HHV4_", "", gene_id))

ebv_expressed <- quant_ebv_df %>%
    group_by(gene_id, tmp_id) %>%
    filter(mean(TPM > 0) >= 0.5) %>%
    ungroup()
   
ebv_gids <- ebv_expressed %>%
    distinct(gene_id) %>%
    mutate(gene_id = sub("HHV4_", "", gene_id),
	   gid = gene_id,
	   gid = case_when(grepl("B*RF1", gene_id) ~ sub("\\.[12]$", "", gene_id),
			   gene_id == "EBNA-1.2" ~ sub("\\.2$", "", gene_id),
			   gene_id == "EBNA-3B/EBNA-3C" ~ sub("/", ":", gene_id),
			   TRUE ~ gene_id)) %>%
    separate_rows(gid, sep = ":")

ebv_bed <- ebv_expressed %>%
    left_join(ebv_gids, by = "gene_id") %>%
    left_join(ebv_annots, by = c("gid" = "gene_id")) %>%
    unite(gene_id, c("gene_id", "tmp_id"), sep = "_") %>%
    select(chr, start, end, id = gene_id, gid, strd = strand, sampleid, tpm = TPM) %>%
    pivot_wider(names_from = sampleid, values_from = tpm) %>%
    select(names(bed_humans))
    
bed_final <- bind_rows(bed_humans, ebv_bed) %>%
    distinct(chr, id, .keep_all = TRUE) %>%
    rename(`#chr` = chr)

write_tsv(bed_final, "./geuvadis_salmon_quants_ebv.bed")


geuvadis_ids <- 
    "https://ftp.ebi.ac.uk/biostudies/fire/E-GEUV-/001/E-GEUV-1/Files/E-GEUV-1.sdrf.txt" |>
    read_tsv() |>
    distinct(id_kgp = `Source Name`, 
             id_geuv = `Comment[ENA_RUN]`,
             pop = `Characteristics[ancestry category]`)

quant_gene_df |>
    left_join(gene_annots, join_by(gene_id)) |>
    filter(gene_name == "IL12A") |>
    left_join(geuvadis_ids, join_by(sampleid == id_geuv)) |>
    arrange(id_kgp) |>
    pull(id_kgp)

### SNRPC
snrpc <- 
    gene_annots |> 
    filter(gene_name == "SNRPC") |>
    pull(gene_id)

snrpc_txs <- 
    transcript_annots |>
    filter(gene_id == snrpc) |>
    pull(transcript_id)

quant_df <- 
    list.files("./results/salmon", pattern = "quant\\.sf", full.names=TRUE, recursive = TRUE) |>
    {function(x) setNames(x, str_extract(x, "ERR\\d+"))}() |>
    future_map_dfr(~read_tsv(.) |> filter(Name %in% snrpc_txs), .id = "sampleid")

