library(tidyverse)


# Save IDs of Europeans in 1000G
dat <- "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index" |>
    read_tsv(comment = "##")

europeans <- c("CEU", "TSI", "GBR", "FIN", "IBS")

dat |>
    distinct(SAMPLE_NAME, POPULATION) |>
    filter(POPULATION %in% europeans) |>
    arrange(SAMPLE_NAME) |>
    pull(SAMPLE_NAME) |>
    write_lines("./data/europeans_samples.txt")


# Save regions in Bentham et al
regions <- read_tsv("../data/coloc_input/regions_bentham.tsv") |>
    mutate(coord = paste0("chr", coord))

write_lines(regions$coord, "./data/regions_bentham.txt")


# susie test
library(susieR)

data("N3finemapping")

Rin <- cor(N3finemapping$X)

dim(Rin)

Rin[1:10, 1:10]
