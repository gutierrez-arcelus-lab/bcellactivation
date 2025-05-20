library(tidyverse)
library(readxl)

dat <- 
    "./common_data/3. dna_antibody_nodate-deidentified.xlsx" |>
    read_excel() |>
    janitor::clean_names() |>
    select(subject_id, age_test, group_id, test_id, result, 
	   comment = result_text, ref_units = reference_units, ref_range = reference_range)

dat |> count(test_id, sort = TRUE) |> print(n = Inf)

test_1 <- dat |>
    filter(test_id %in% c("1-303", "DNAAB", "5200004438"),
	   ref_units %in% c("IU", "IU/mL")) |>
    separate(ref_range, c("ref_min", "ref_thres"), sep = "-", convert = TRUE) |>
    mutate(result = sub("^>|<", "", result),
	   result = as.numeric(result)) |>
    mutate(positive = as.integer(result > ref_thres)) |>
    select(subject_id, positive)

test_2 <- dat |> 
    filter(test_id %in% c("SQ-DNA", "5200004436")) |>
    mutate(positive = case_when(grepl("Positive", result, ignore.case = TRUE) ~ 1L,
			      grepl("Negative", result, ignore.case = TRUE) ~ 0L,
			      .default = NA_integer_)) |>
    select(subject_id, positive)
    
# Not valid results
dat |> 
    filter(test_id == "1-304")

test_3 <- dat |>
    filter(test_id == "DSDNA1", 
	  ref_units %in% c("IU", "IU/mL")) |>
    mutate(result = sub("^<|>", "", result),
	   result = as.numeric(result),
	   positive = as.integer(result > 75)) |>
    select(subject_id, positive)

test_4 <- dat |> 
    filter(test_id == "DNA") |>
    mutate(positive = case_when(result == "POS" ~ 1L,
				result == "NEG" ~ 0L,
			      .default = NA_integer_)) |>
    select(subject_id, positive)
    
valid_tests <- 
    bind_rows(test_1, test_2, test_3, test_4) |>
    drop_na(positive) |>
    group_by(subject_id) |>
    summarise(ever_positive = as.integer(any(positive == 1L))) |>
    ungroup()

write_tsv(valid_tests, "./mgb_data/results_ab.tsv")
