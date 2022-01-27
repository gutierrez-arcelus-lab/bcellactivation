library(data.table)
library(tidyverse)
library(qvalue)

annotations <- "../data/gencode.v38.primary_assembly.annotation.gtf" %>%
    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc") %>%
    filter(X3 == "exon") %>%
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) %>%
    select(chr = 1, start = 4, end = 5, gene_id, gene_name)

ase <- read_tsv("./results/gatk/ERR188022.ase.counts") %>%
    mutate(p = map2_dbl(refCount, totalCount, ~binom.test(.x, .y, p = .5, alternative = "two.sided")$p.value)) %>%
    select(chr = 1, pos = 2, ref = 4, alt = 5, ref_n = 6, total = 8, p) %>%
    mutate(q = qvalue(p)$qvalues)

ase_wasp <- read_tsv("./results/gatk/ERR188022.ase.wasp.counts") %>%
    mutate(p = map2_dbl(refCount, totalCount, ~binom.test(.x, .y, p = .5, alternative = "two.sided")$p.value)) %>%
    select(chr = 1, pos = 2, ref = 4, alt = 5, ref_n = 6, total = 8, p) %>%
    mutate(q = qvalue(p)$qvalues)

ase_wasp_dt <- ase_wasp %>%
    mutate(pos2 = pos) %>%
    select(chr, pos, pos2) %>%
    as.data.table()

annots_dt <- annotations %>%
    select(chr, start, end) %>%
    as.data.table()

setkey(annots_dt)

ase_wasp_annot <- foverlaps(ase_wasp_dt, annots_dt, 
	  by.x = names(ase_wasp_dt),
	  type = "within", 
	  mult = "all")

annotations %>%
    arrange(chr, start) %>%
    filter(chr == "chr1", start >= 14000)

write_tsv(ase, "./results/gatk/ERR188022.ase.pvalues.tsv")
write_tsv(ase_wasp, "./results/gatk/ERR188022.ase.wasp.pvalues.tsv")

ase_qtltools <- read_tsv("./results/qtltools/ERR188022.ase") %>%
    select(chr = CHR, pos = POS, p = PVALUE) %>%
    mutate(q = qvalue(p)$qvalues)

write_tsv(ase_qtltools, "./results/qtltools/ERR188022.ase.fdr.txt")

ase_df <- 
    bind_rows("asereadcounter" = select(ase, chr, pos, p, q),
	      "asereadcounter_wasp" = select(ase_wasp, chr, pos, p, q),
	      "qtltools" = ase_qtltools, 
	      .id = "method")

write_tsv(ase_df, "./ase_data.tsv")


bind_rows("all" = count(ase_df, method),
	  "FDR = 5%" = ase_df %>% filter(q <= 0.05) %>% count(method), 
	  .id = "category") %>%
write_tsv("./plots_eda/summ_total_asevars.tsv")
    
