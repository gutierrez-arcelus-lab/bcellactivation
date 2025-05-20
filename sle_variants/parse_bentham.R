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

summ_stats_original |>
    filter(chrom == 3, pos > 159e6, pos < 161e6) |>
    arrange(p) |>
    select(chrom:beta)


summ_stats_original |> filter(rsid %in% balloons1) |> arrange(p)
summ_stats_original |> filter(rsid %in% balloons2) |> arrange(p)

summ_stats |>
    filter(chromosome == 17, base_pair_location > 35056109 - 500, base_pair_location < 35056109 + 500)

summ_stats |>
    filter(chromosome == 14, base_pair_location == 56820049)


# WANG
wang <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Wang2021/ASN/GCST90011866_buildGRCh37.tsv" |>
    data.table::fread() |>
    as_tibble()

wang |> filter(variant_id == "rs564799")


# Compare versions of Bentham et al
gwascat_harm <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/harmonized/26502338-GCST003156-EFO_0002690-build37.f.tsv.gz" |>
    data.table::fread() |>
    as_tibble() |>
    select(chrom = chromosome, pos = base_pair_location, rsid = variant_id, other_allele, effect_allele, beta)

open_gwas <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/opengwas/ebi-a-GCST003156.vcf.gz" |>
    data.table::fread(skip = "#CHROM") |>
    as_tibble() |>
    select(chrom = `#CHROM`, pos = POS, rsid = ID, REF, ALT, FORMAT, stats = starts_with("EBI")) 

open_gwas_stats <- 
    open_gwas |>
    select(FORMAT, stats) |>
    separate_rows(c(FORMAT, stats), sep = ":") |>
    pivot_wider(names_from = FORMAT, values_from = stats)


# Compare with Wang
bentham <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Bentham/harmonized/26502338-GCST003156-EFO_0002690-build37.f.tsv.gz" |>
    data.table::fread() |>
    as_tibble() |>
    mutate(z = beta/standard_error) |>
    select(chrom = chromosome, pos = base_pair_location, rsid = variant_id, other_allele, effect_allele, z)

wang <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Wang2021/EUR/Meta_Results.txt" |>
    data.table::fread() |>
    as_tibble() |>
    select(chrom = CHR, pos = BP, rsid = SNP, other_allele = A1lele1, effect_allele = Allele2, z = Zscore) |>
    mutate_at(vars(ends_with("allele")), toupper) |>
    mutate(chrom = as.character(chrom))

bentham_wang_df <- 
    inner_join(bentham, wang, join_by(chrom, pos, rsid), suffix = c("_bentham", "_wang")) |>
    mutate(chrom = paste0("chr", chrom),
	   chrom = factor(chrom, levels = paste0("chr", 1:22))) |>
    arrange(chrom, pos)

bentham_plot <- 
    ggplot(bentham_wang_df |> sample_frac(.1), 
	   aes(x = abs(z_bentham), y = abs(z_wang))) +
    geom_abline() +
    geom_point(size = .25, alpha = .25) +
    facet_wrap(~chrom, nrow = 5) +
    theme_bw() +
    theme(panel.grid = element_blank())

ggsave("./bentham_wang.png", bentham_plot)


