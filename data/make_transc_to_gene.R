library(tidyverse)

annotations <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens",
	      "gencode.v38.primary_assembly.annotation.gtf") %>%
    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc") %>%
    filter(X3 == "transcript", X1 %in% paste0("chr", c(1:22, "X"))) %>%
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
           gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	   transcript_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"),
           transcript_type = str_extract(X9, "(?<=transcript_type\\s\")[^\"]+")) %>%
    select(gene_id, gene_name, transcript_id, transcript_type)

write_tsv(annotations, "./transcript_annots_v38.tsv")
