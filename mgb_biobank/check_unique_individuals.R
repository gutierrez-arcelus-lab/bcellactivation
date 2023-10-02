library(tidyverse)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 5)

mgbdir <- "/lab-share/IM-Gutierrez-e2/Public/External_datasets/mgb_biobank"  

batches <- list.files(mgbdir, pattern = "^04")

ids_df <- 
    file.path(mgbdir, batches, "vcf/chr22.dose.vcf.gz") %>%
    map_df(. %>% 
	   read_tsv(comment = "##", n_max = 1) %>%
	   select(-(1:9)) %>%
	   names() %>%
	   tibble(id = .))

ids_df %>%
    separate(id, c("dummy", "id"), sep = "-") %>%
    count(id, sort = TRUE)

# test
donors_tae <- readxl::read_excel("/home/ch229163/24 donrs ID with genotype for rs117701653.xlsx") %>%
    mutate_at(vars(Subject_ID), as.character)


ids_df <- ids_df %>%
    separate(id, c("dummy", "id"), sep = "-")

all(donors_tae$Subject_ID %in% ids_df$id)
