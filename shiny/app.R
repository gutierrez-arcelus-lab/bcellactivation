library(shiny)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)

## Function
format_timecourse <- function(expression_data, chosen_gene) {

    d <- 
        expression_data |> 
        dplyr::filter(gene_label == chosen_gene) |>
        tidyr::unite("condition", c(stim, timep), sep = "_", remove = FALSE)

    d_tmp <- dplyr::filter(d, condition != "Unstim_0")

    dplyr::filter(d, condition == "Unstim_0") |> 
        dplyr::select(-stim) |>
        tidyr::expand_grid(dplyr::distinct(d_tmp, stim)) |>
        dplyr::bind_rows(d_tmp)
}  

## Metadata
plot_colors <- 
    readr::read_tsv("./data/figure_colors_v2.txt", 
                    col_names = c("stim", "timep", "col")) |>
    tidyr::unite("lab", c(stim, timep), sep = "_") |>
    tibble::deframe()

atac_colors <- plot_colors[c("Unstim_0", "Unstim_24", "IL-4c_24", "TLR7c_24", "BCRc_24", "DN2c_24")]
atac_colors[c("TLR7c_24", "BCRc_24", "DN2c_24")] <- plot_colors[c("TLR7c_24", "BCRc_48", "DN2c_48")] 

## Timecourse data
gene_list <- readr::read_lines("./data/genes.txt")

gene_exp <- readr::read_rds("./data/deseq_normalized_counts.rds")

## Gene tracks for ATAC-seq plot
ah <- AnnotationHub::AnnotationHub()
ens_data <- ah[["AH98047"]]

## Atac-seq tracks
bigwigs <- 
    "./data/bigwig" |>
    list.files(full.name = TRUE, pattern = "*.bigWig")

names(bigwigs) <- sub("^([^_]+_\\d+).+$", "\\1", basename(bigwigs))

## Cite-seq data
sc_data <- SeuratDisk::LoadH5Seurat("./data/bcells.h5Seurat")

sc_meta <- 
    sc_data@meta.data |>
    as.data.frame() |>
    tibble::rownames_to_column("barcode") |>
    tibble::as_tibble() |>
    dplyr::select(barcode, hto = dmm_hto_call)

umap_df <-
    Seurat::Embeddings(sc_data, reduction = "umap") |>
    as.data.frame() |>
    tibble::rownames_to_column("barcode") |>
    tibble::as_tibble() |>
    dplyr::left_join(sc_meta, dplyr::join_by(barcode)) |>
    dplyr::mutate(hto = dplyr::recode(hto, 
                                      "Unstim 0h" = "Unstim_0",
                                      "IL4 24h" = "IL-4c_24",
                                      "IL4 72h" = "IL-4c_72",
                                      "TLR7 24h" = "TLR7c_24",
                                      "TLR7 72h" = "TLR7c_72",
                                      "BCR 24h" = "BCRc_24",
                                      "BCR 72h" = "BCRc_72",
                                      "DN2 24h" = "DN2c_24",
                                      "DN2 72h" = "DN2c_72"))

umap_stims <-
    ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(fill = hto), size = 2, 
               shape = 21, stroke = .05, color = "black") +
    scale_fill_manual(values = plot_colors) +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          panel.grid = element_blank()) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    guides(fill = guide_legend(title = "Stim:",
                               title.position = "top",
                               override.aes = list(size = 3)))

## App
ui <- navbarPage(
    "B cell activation",
    tabPanel("Bulk RNA-seq",
             selectizeInput("generna", "Select gene:", choices = NULL),
             plotOutput("plotrna", width = "100%", height = "200px")
             ),
    tabPanel("Bulk ATAC-seq",
             selectizeInput("geneatac", "Select gene:", choices = NULL),
             numericInput("atacwindow", "Window Size", value = 10e3, min = 1e3, max = 500e3),
             plotOutput("plotatac", width = "100%", height = "600px")
             ),
    tabPanel("Single-cell RNA-seq",
             selectizeInput("genesc", "Select gene:", choices = NULL),
             plotOutput("plotsc", width = "100%", height = "600px")
    )
)

server <- function(input, output, session) {

    ## Bulk RNA
    updateSelectizeInput(session, 
                         "generna", 
                         choices = gene_list, 
                         server = TRUE)
    
    datarna <- reactive(gene_exp |> format_timecourse(input$generna))
    
    output$plotrna <- 
        renderPlot(ggplot(data = datarna()) +
                       geom_quasirandom(aes(x = timep, y = norm_counts, fill = condition),
                                        size = 2.5, stroke = .25, shape = 21,
                                        method = "smiley", width = .2) +
                       scale_fill_manual(values = plot_colors) +
                       facet_grid(cols = vars(stim), space = "free", scale = "free_x") +
                       theme_bw() +
                       guides(fill = "none") +
                       labs(x = "hours", y = "Norm. counts"), 
                   res = 96)
    
    ## Bulk ATAC
    updateSelectizeInput(session, 
                         "geneatac", 
                         choices = gene_list, 
                         server = TRUE)
    loc <-
        reactive(
            locuszoomr::locus(gene = input$geneatac, 
                              flank = input$atacwindow,
                              ens_db = ens_data)
        )
    
    atac_peaks <-
        reactive({

            interv <- 
                GenomicRanges::GRanges(paste0("chr", loc()$seqname), 
                                       IRanges(loc()$xrange[1], loc()$xrange[2]))
            
            atac_ranges <- purrr::map(bigwigs, ~rtracklayer::import(., which = interv))
            
            atac_covered <-
                atac_ranges |>
                purrr::map_dfr(as.data.frame, .id = "stim") |>
                dplyr::select(stim, start, end, score)
            
            atac_gaps <- 
                atac_ranges |>
                purrr::map_dfr(~ranges(.) |> gaps() |> as.data.frame(), .id = "stim") |>
                dplyr::mutate(score = 0) |>
                dplyr::select(stim, start, end, score)
                
            dplyr::bind_rows(atac_covered, atac_gaps) |>
                dplyr::mutate(stim = stringr::str_replace(stim, "unst", "Unstim"),
                              stim = stringr::str_replace(stim, "IL4", "IL-4c"),
                              stim = stringr::str_replace(stim, "TLR7", "TLR7c"),
                              stim = stringr::str_replace(stim, "BCR", "BCRc"),
                              stim = stringr::str_replace(stim, "DN2", "DN2c"),
                              stim = factor(stim, levels = names(atac_colors))) |>
                dplyr::arrange(stim, start) |>
                tidyr::pivot_longer(start:end, names_to = "dummy", values_to = "pos")
        })
    
    plot_atac <-
        reactive({

            ggplot(atac_peaks()) +
            geom_ribbon(aes(x = pos, ymin = 0, ymax = score, color = stim, fill = stim),
                        linewidth = .5, outline.type = "full", alpha = .5) +
            scale_x_continuous(limits = loc()$xrange, 
                               labels = function(x) round(x/1e6L, 2),
                               expand = c(0, 0)) +
            scale_color_manual(values = atac_colors) +
            scale_fill_manual(values = atac_colors) +
            facet_wrap(~stim, ncol = 1, strip.position = "right") +
            theme_minimal() +
            theme(
                axis.text = element_text(size = 12),
                legend.position = "none",
                strip.text.y.right = element_text(angle = 0, size = 12),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                plot.background = element_rect(color = "white", fill = "white")) +
            labs(x = NULL)
        })
        
    gene_tracks <- 
        reactive({
            locuszoomr::gg_genetracks(loc(), cex.text = 1) + 
            scale_x_continuous(limits = loc()$xrange/1e6,
                               labels = function(x) round(x, 2),
                               expand = c(0, 0)) +
            theme_minimal() +
            theme(
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 12),
                plot.background = element_rect(color = "white", fill = "white"))
        })
    
    output$plotatac <- renderPlot(plot_atac() / gene_tracks() + plot_layout(heights = c(1, .3)))
    
    #Single-cell RNA-seq
    updateSelectizeInput(session, 
                         "genesc", 
                         choices = gene_list, 
                         server = TRUE)
    
    sc_data_gene <- 
        reactive({
            Seurat::GetAssayData(sc_data, slot = "data")[input$genesc, , drop=FALSE] |> 
                t() |>
                as.data.frame() |> 
                tibble::rownames_to_column("barcode") |>
                tibble::as_tibble() |> 
                dplyr::select(barcode, value = 2) |>
                dplyr::right_join(umap_df, dplyr::join_by(barcode)) |>
                dplyr::arrange(value) |> 
                dplyr::mutate(barcode = forcats::fct_inorder(barcode))
            })
    
    umaps_plot_expr <-
        reactive(
            ggplot(data = sc_data_gene(), 
                   aes(x = UMAP_1, y = UMAP_2)) +
                geom_point(aes(color = value), size = 2, stroke = 0) +
                scale_color_gradient2(low = "grey70", 
                                      mid = "lightyellow", 
                                      high = "#ff0000",
                                      midpoint = mean(sc_data_gene()$value[sc_data_gene()$value > 0])) +
                theme_minimal() +
                theme(axis.text = element_text(size = 12),
                      axis.title = element_text(size = 12),
                      legend.text = element_text(size = 12),
                      legend.title = element_text(size = 12),
                      panel.grid = element_blank()) +
                guides(color = guide_colorbar(barwidth = .2, barheight = 8))
        )
    
    output$plotsc <- 
        renderPlot(umap_stims + 
                       plot_spacer() + 
                       umaps_plot_expr() + 
                       plot_layout(widths = c(1, .1, 1)))
}

shinyApp(ui, server)