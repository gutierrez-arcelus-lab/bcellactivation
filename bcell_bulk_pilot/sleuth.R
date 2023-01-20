library(tidyverse)
library(sleuth)

sample_ids <- read_lines("./bcell_samples.txt")
results_path <- file.path("./bcell_quant/kallisto", sample_ids) 
conditions <- str_extract(sample_ids, "(\\d{2}hr_[^_]+)")

info_df <- tibble(sample = sample_ids, condition = conditions, path = results_path)  

transc_to_gene <- read_tsv("./data/transc_to_gene.tsv")

# construct the sleuth object
sleuth_obj <- sleuth_prep(info_df, target_mapping = transc_to_gene)

# fit model
sleuth_obj <- sleuth_fit(sleuth_obj, ~condition, 'full')
sleuth_obj <- sleuth_fit(sleuth_obj, ~1, 'reduced')
sleuth_obj <- sleuth_lrt(sleuth_obj, 'reduce', 'full')

#so <- sleuth_wt(sleuth_obj, which_beta = "72hr_IgG") 

