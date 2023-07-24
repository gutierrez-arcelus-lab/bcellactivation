library(tidyverse)
library(readxl)

dat <- 
    "./common_data/3. dna_antibody_nodate-deidentified.xlsx" |>
    read_excel() |>
    janitor::clean_names() |>
    select(subject_id, age_test, group_id, test_id, result, 
	   comment = result_text, ref_units = reference_units, ref_range = reference_range)

dat |> count(test_id, sort = TRUE)

dat |> 
    filter(ref_range == "0.0-0.0") |>
    count(test_id, ref_units, sort = TRUE)


