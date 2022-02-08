library(data.table)
library(tidyverse)
library(qvalue)

annotations <- "../data/gencode.v38.primary_assembly.annotation.gtf" %>%
    read_tsv(comment = "#", col_names = FALSE, col_types = "ccciicccc") %>%
    filter(X3 == "exon") %>%
    mutate(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+")) %>%
    select(chr = 1, start = 4, end = 5, gene_id, gene_name)

sampleids <- read_lines("./bcell_samples.txt")


aseread <- "./results/ase/%s.asereadcounter.txt" %>%
    sprintf(sampleids) %>%
    setNames(sampleids) %>%
    map_df(. %>% read_tsv() %>% 
	   mutate(p = map2_dbl(refCount, totalCount, 
			       ~binom.test(.x, .y, p = .5, alternative = "two.sided")$p.value)) %>%
	   select(chr = 1, pos = 2, ref = 4, alt = 5, ref_n = 6, total = 8, p) %>%
	   mutate(q = qvalue(p)$qvalues), 
           .id = "id")

aseread_dt <- aseread %>%
    distinct(chr, pos) %>%
    mutate(pos2 = pos) %>%
    as.data.table()

annots_dt <- annotations %>%
    select(chr, pos = start, pos2 = end) %>%
    as.data.table()

setkey(annots_dt)

annot_positions <- 
    foverlaps(aseread_dt, annots_dt, 
	      by.x = names(aseread_dt),
	      type = "within", 
	      mult = "all") %>%
    as_tibble() %>%
    select(chr, start = pos, end = pos2, pos = i.pos) %>% 
    left_join(annotations, by = c("chr", "start", "end")) %>%
    distinct(chr, pos, gene_id, gene_name) %>%
    rowwise() %>%
    mutate(annot = case_when(!is.na(gene_id) & !is.na(gene_name) ~ paste(gene_id, gene_name, sep = ":"),
			     TRUE ~ NA_character_)) %>%
    group_by(chr, pos) %>%
    summarise(annot = paste(annot, collapse = ";")) %>%
    ungroup()
    
aseread_annotated <- aseread %>%
    left_join(annot_positions, by = c("chr", "pos")) %>%
    add_column(method = "ASEReadCounter", .before = 1)


ase_qtltools <- "./results/ase/%s.subsamp%d.qtltools.ase" %>%
    sprintf(sampleids, rep(c(25, 50, 75), each = 5)) %>%
    setNames(paste0(sampleids, "-QTLtools.", rep(c(25, 50, 75), each = 5))) %>%
    map_df(. %>% read_tsv() %>%
	   select(chr = CHR, pos = POS, ref = REF_ALLELE, alt = ALT_ALLELE,
		  ref_n = REF_COUNT, total = TOTAL_COUNT, p = PVALUE, 
		  annot = EXON_INFO) %>%
	   mutate(q = qvalue(p)$qvalues) %>%
	   extract(annot, c("gene_id", "gene_name"), "(ENSG[^:]+):.+:.+:(.+)") %>%
	   rowwise() %>%
	   mutate(annot = case_when(!is.na(gene_id) & !is.na(gene_name) ~ paste(gene_id, gene_name, sep = ":"),
				    TRUE ~ NA_character_)) %>%
	   select(chr, pos, ref, alt, ref_n, total, p, q, annot),
           .id = "id") %>%
    separate(id, c("id", "method"), sep = "-") %>%
    select(method, id, everything())

out <- bind_rows(aseread_annotated, ase_qtltools) 

write_tsv(out, "./results/ase/ase_compiled.tsv")

