library(tidyverse)

regions <- 
    "./data/pt_5e-4_regions.tsv" |>
    read_tsv(col_names = c("rsid", "region_coord")) |>
    rowid_to_column("region")

tag_files <- 
    system("echo $TEMP_WORK", intern = TRUE) |>
    file.path("vcf/prs/hg19/pt_5e-4/%s_tagvariants.tsv") |>
    sprintf(regions$rsid) |>
    setNames(regions$region)
    
out <- tag_files |>
    map(~read_tsv(., col_types = "iiccccd")) |>
    keep(~nrow(.) > 0) |>
    bind_rows(.id = "region") |>
    mutate(region = as.integer(region),
	   var_type = factor(var_type, levels = c("lead", "tag"))) |>
    arrange(region, var_type, chr, pos, r2)

# Three sentinel variants are missing because they are not polymorphic in the 1000 Genomes subset
# Add them from original summ stats
sentinels <- 
    "/lab-share/IM-Gutierrez-e2/Public/GWAS/SLE/Khunsriraksakul/Khunsriraksakul_summstats_mamt.txt" |>
    read_tsv() |>
    select(rsid, chr, pos = `pos(hg19)`, ref = `effect allele`, alt = `other allele`)

missing_df <- anti_join(regions, out, join_by("region")) |>
    left_join(sentinels) |>
    mutate(var_type = "lead", r2 = 1) |>
    select(region, chr, pos, rsid, var_type, ref, alt, r2)

bind_rows(out, missing_df) |>
    arrange(region, var_type, chr, pos, r2) |>
    write_tsv("./data/pt_5e-4_lead_and_tag.tsv")


