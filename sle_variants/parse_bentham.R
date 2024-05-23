library(tidyverse)

bentham_regions <-
    "paper_data/bentham_tab1.tsv" |>
    read_tsv() |>
    mutate(start = pos - 5e5L, end = pos + 5e5L) |>
    select(chr, locus, start, end)

summ_stats <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/harmonized/26502338-GCST003156-EFO_0002690-build37.f.tsv.gz" |>
    data.table::fread() |>
    as_tibble()

summ_stats_regions <-
    inner_join(summ_stats, 
	       bentham_regions, 
	       join_by(chromosome == chr, between(base_pair_location, start, end)))

balloons1 <- 
    c("rs7285053", "rs8078864", "rs17091347", "rs73050535", "rs1034009", "rs10264693", 
      "rs17168663", "rs16895550", "rs9969061", "rs11962557", "rs13170409", "rs6532924", 
      "rs17087866", "rs12081621", "rs4074976", "rs4661543", "rs13019891", "rs512681", 
      "rs10200680", "rs2573219","rs9852014", "rs1464446", "rs1078324", "rs7386188", 
      "rs7823055", "rs11928304", "rs12309414", "rs12948819")
 
balloons2 = 
    c("rs34703115","rs55684314","rs77601723","rs73068668","20:7357645:i","rs111478576",
      "rs1120809","rs11905662","rs11907821","rs2423194","rs34995148","rs56412650",
      "rs57830939","rs58072943","rs58605847","rs58829623","rs60461307","rs7262084",
      "rs73894580","rs73894583","rs73894585","rs73894599","rs73897155","rs77222522")

length(balloons1)
length(balloons2)

summ_stats_regions |>
    filter(variant_id %in% balloons2)

summ_stats |> filter(variant_id %in% balloons1)
summ_stats |> filter(variant_id %in% balloons2)

summ_stats_original <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz" |>
    data.table::fread() |>
    as_tibble()

summ_stats_original |> filter(rsid %in% balloons1) |> arrange(p)
summ_stats_original |> filter(rsid %in% balloons2) |> arrange(p)

summ_stats |>
    filter(chromosome == 17, base_pair_location > 35056109 - 500, base_pair_location < 35056109 + 500)

summ_stats |>
    filter(chromosome == 14, base_pair_location == 56820049)
