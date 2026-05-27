library(tidyverse)
library(glue)
library(cowplot)
library(patchwork)

my_stims <- c("CD40L", "TLR9", "TLR7", "BCR", "BCR-TLR7", "DN2")
 
go_all <-
    glue("../bcell_lowinput/wgcna/data/{my_stims}_go.tsv") |>
    setNames(my_stims) |>
    map_dfr(read_tsv, .id = "stim") |>
    mutate(log_p = -log10(pvalue),
           gene_r = map_dbl(GeneRatio, ~eval(parse(text = .)))) |>
    select(stim, module, Description, gene_r, log_p)
 
dummy_plot <-
    ggplot(go_all, aes(x = 1, y = 1)) +
    geom_point(aes(fill = log_p, size = gene_r), shape = 21, stroke = .2) +
    scale_size(limits = range(go_all$gene_r),
               range = c(.5, 2.5)) +
    scale_fill_gradient(low = "white", high = "black",
                        limits = c(0, max(go_all$log_p))) +
    guides(
        fill = guide_colorbar("logP:", barheight = 5, barwidth = .25),
        size = guide_legend("Gene\nRatio:")
    ) +
    theme_minimal() +
    theme(
        legend.position = "right",
        legend.box = "vertical",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.margin = margin(l = 10, unit = "pt"),
    )

global_legend <- get_legend(dummy_plot)

build_module_row <- function(mod_name, modules_df, kim_df, go_res) {
    
    mod_data <- filter(modules_df, module_ix == mod_name)
    kim_data <- filter(kim_df, module_ix == mod_name)
    go_data  <- 
        filter(go_res, module == mod_name) |>
        mutate(Description = fct_inorder(as.character(Description)))
        
    
    p1 <- 
        ggplot(mod_data, aes(x = time, y = value)) + 
        geom_line(aes(group = module_ix), linewidth = .75) +
        scale_y_continuous(limits = range(modules_df$value),
                           breaks = range(modules_df$value)) +
        facet_wrap(~module_ix, strip.position = "left") +
        theme_minimal() +
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.length.x = unit(0, "pt"),
              panel.grid.major.x = element_line(linewidth = .2),
              panel.grid.major.y = element_line(linewidth = .2),
              panel.grid.minor.y = element_blank(),
              strip.placement = "outside",
              strip.text.y.left = element_text(angle = 0, size = 5, margin = margin(l = 0, r = 0)),
              strip.clip = "off",
              plot.margin = margin(t = 0, b = -1, unit = "lines"),
        )
    
    p2 <- 
        ggplot(kim_data) +
        geom_point(aes(x = " ",
                       fct_rev(fct_inorder(str_pad(gene_name, width = 15, side = "right"))),
                       fill = kim), 
                   shape = 21, size = 2, stroke = .1,
                   show.legend = FALSE) +
        scale_y_discrete(position = "right") +
        scale_fill_continuous(low = "beige", high = "firebrick",
                              limits = range(kim_df$kim)) +
        theme_minimal() + 
        theme(axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y.right = element_text(size = 4, hjust = 0, face = "italic", family = "mono"),
              axis.ticks.x = element_blank(),
              axis.ticks.length.x = unit(0, "pt"),
              panel.grid = element_blank(),
              plot.margin = margin(t = 0, b = -1, unit = "lines"),
        ) +
        coord_cartesian(clip = "off")
    
    p3 <- 
        ggplot(go_data, aes(x = " ", y = Description)) +
        geom_point(aes(fill = -log10(pvalue), size = gene_r), 
                   shape = 21, stroke = .1, show.legend = FALSE) +
        scale_y_discrete(position = "right", limits = rev) +
        scale_size(limits = range(go_all$gene_r),
                   range = c(.5, 2.5)) +
        scale_fill_gradient(low = "white", high = "black",
                            limits = c(0, max(go_all$log_p))) +
        theme_minimal() +
        theme(
              axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 4, family = "mono"),
              axis.ticks.x = element_blank(),
              axis.ticks.length.x = unit(0, "pt"),
              panel.grid = element_blank(),
              plot.margin = margin(t = 0, b = -1, unit = "lines"),
        ) +
        coord_cartesian(clip = "off")
    
    # 5. Combine into a single, perfectly aligned row
    row_plot <- 
        p1 + p2 + p3 + 
        plot_layout(widths = c(2, 0.1, .2)) &
        theme(plot.margin = margin(t = 0, r = .25, b = 1, l = .25, unit = "pt"))
    
    return(row_plot)
}

build_stim_panel <- function(stim_i) {
    
    module_sizes <-
        glue("../bcell_lowinput/wgcna/data/{stim_i}_modules.tsv") |>
        read_tsv() |>
        count(module) |>
        arrange(desc(n)) |>
        filter(module != "grey") |>
        rowid_to_column("ix") |>
        select(module, ix) |>
        mutate(ix = paste("Module", ix)) |>
        mutate(ix = fct_inorder(ix))
    
    eigengenes_df <-
        glue("../bcell_lowinput/wgcna/data/{stim_i}_eigen.tsv") |>
        read_tsv() |>
        select(-MEgrey) |>
        separate(sample_name, 
    	     c("donor_id", "stim", "time"), 
    	     sep = "_") |>
        pivot_longer(starts_with("ME"), names_to = "module") |>
        mutate(module = str_remove(module, "^ME"),
    	   time = parse_number(time),
    	   time = factor(time, levels = sort(unique(time)))) |>
        left_join(module_sizes, join_by(module)) |>
        select(donor_id, module_ix = ix, module, time, value) |>
        arrange(module_ix, donor_id, time)
    
    modules_df <-
        eigengenes_df |>
        group_by(module, module_ix, time) |>
        summarise(value = mean(value)) |>
        ungroup() |>
        arrange(module_ix, time)
    
    kim_df <- 
        glue("../bcell_lowinput/wgcna/data/{stim_i}_kim.tsv") |>
        read_tsv() |>
        left_join(module_sizes, join_by(module)) |> 
        select(gene_id, gene_name, module, module_ix = ix, kim) |>
        arrange(module_ix, desc(kim)) |>
        group_by(module) |>
        top_n(5, kim) |>
        ungroup()
    
    go_res <-
        glue("../bcell_lowinput/wgcna/data/{stim_i}_go.tsv") |>
        read_tsv() |>
        group_by(module) |>
        slice_max(-log10(pvalue), n = 5, with_ties = FALSE) |>
        ungroup() |>
        mutate(Description = str_trunc(Description, width = 36),
               gene_r = map_dbl(GeneRatio, ~eval(parse(text = .)))) |>
        right_join(module_sizes, join_by(module)) |>
        select(module = ix, Description, pvalue, gene_r) |>
        arrange(module, pvalue) |>
        mutate(Description = fct_inorder(Description))
        
    modules <- unique(module_sizes$ix)
    
    all_rows <- lapply(modules, function(mod) {
        build_module_row(mod, modules_df, kim_df, go_res)
    })
    
    stim_figure <- wrap_plots(all_rows, ncol = 1)
    
    wrap_elements(panel = stim_figure) +
        ggtitle(recode(stim_i, "CD40L" = "CD40c", "TLR9" = "TLR9c", "TLR7" = "TLR7c", "BCR" = "BCRc", "BCR-TLR7" = "BCR/TLR7c", "DN2" = "DN2c")) +
        theme(plot.title = element_text(face = "bold", size = 6, hjust = .5, margin = margin(b = 0, unit = "pt")),
              plot.margin = margin(t = 1, r = .25, l = .25, b = 0, unit = "pt")
        ) +
        coord_cartesian(clip = "off")
}

all_stim_panels <- lapply(my_stims, build_stim_panel)

final_figure <- 
    wrap_plots(all_stim_panels, ncol = 2, heights = c(1, 1, 1.2))

final_out <- 
    plot_grid(final_figure, global_legend, nrow = 1, rel_widths = c(1, 0.1)) +
    theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("./modules_all_stims.png", final_out, width = 6.5, height = 8, dpi = 300)