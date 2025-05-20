library(tidyverse)
library(ggforce)
library(VennDiagram)
library(RColorBrewer)

cov <- tibble(depth = as.numeric(read_lines("./depth.txt"))) %>%
    count(depth) %>%
    arrange(depth) %>%
    mutate(p = n/sum(n)) %>%
    mutate(csum = cumsum(p))

ggplot(cov, aes(depth, p)) +
    geom_col() +
    scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1L)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    facet_zoom(depth < 25) +
    labs(y = NULL)

ggsave("./depth.png", height = 3, width = 5)

ggplot(cov, aes(depth, csum)) +
    geom_line() +
    geom_area(aes(fill = "1"), fill = "black", alpha = .5, show.legend = FALSE) +
    scale_y_continuous(labels = function(x) scales::percent(x, accuracy = 1L)) +
    coord_cartesian(xlim = c(0, 15)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    labs(y = "Cummulative % depth")

ggsave("./cumm_depth.png", height = 3, width = 5)

#ASE

ase_df <- read_tsv("./ase_data.tsv")

ase_list <- ase_df %>%
    unite("var", c("chr", "pos"), sep = "_") %>%
    split(.$method) %>%
    map(~pull(., var)) 

mycol <- brewer.pal(3, "Pastel2")

venn.diagram(ase_list, 
             catergory.names = names(ase_list),
             filename = "asevenn.png",
             output = TRUE,
             imagetype = "png",
             fill = mycol,
             lty = "blank")


col_info <- tibble(chr = c(1:22, "X"),
                   i = rep(1:2, 20)[1:23])

ase_formatted <- ase_df %>%
    filter(chr %in% paste0("chr", c(1:22, "X"))) %>%
    mutate(chr = sub("chr", "", chr)) %>%
    left_join(col_info) %>%
    mutate(chr = factor(chr, levels = c(1:22, "X"))) %>%
    arrange(method, chr) %>%
    group_by(method) %>%
    mutate(ix = 1:n()) %>%
    ungroup() 


ggmanhat <- function(data_df, m) {

    plot_data <- data_df %>%
        filter(method == m)

    signif_line <- plot_data %>%
        filter(q <= 0.01) %>%
        arrange(-q) %>%
        slice(1) %>%
        pull(p) %>%
        log10() %>%
        `*`(-1)
    
    ggplot(plot_data, aes(ix, -log10(q))) +
    geom_point(aes(color = factor(i)), 
               alpha = .5,
               size = .5,
               show.legend = FALSE) +
    geom_hline(yintercept = signif_line, linetype = 2) +
    ggsci::scale_color_npg() +
    facet_grid(~chr, 
               scales = "free_x", 
               space = "free_x", 
               switch = "x") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.spacing.x = unit(0, "lines"),
          strip.text = element_text(size = 7)) +
        labs(title = m)
}
 

ggmanhat(ase_formatted, "asereadcounter")
ggsave("./asereadcounter_manhat.png", width = 8)

ggmanhat(ase_formatted, "asereadcounter_wasp")
ggsave("./wasp_manhat.png", width = 8)
   
ggmanhat(ase_formatted, "qtltools")
ggsave("./qtltools_manhat.png", width = 8)







