library(shiny)
library(ggplot2)
library(ggbeeswarm)
library(hdf5r)
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

get_sc_data <- function(gene_i) {
    
    idx <- which(sc_genes == gene_i)
    h5file <- H5File$new("./data/singlecell/bcells_expressed.h5")
    h5data <- h5file[["expr_data"]]
    gene_data <- h5data$read(args = list(idx, quote(expr=)))
    gene_df <- tibble::tibble( "{gene_i}" := gene_data)
    h5file$close_all()
    return(gene_df)
}

## Metadata
plot_colors <- 
    readr::read_tsv("./data/figure_colors_v2.txt", 
                    col_names = c("stim", "timep", "col")) |>
    tidyr::unite("lab", c(stim, timep), sep = "_") |>
    tibble::deframe()

atac_colors <- plot_colors[c("Unstim_0", "Unstim_24", "IL-4c_24", "TLR7c_24", "BCRc_24", "DN2c_24")]
atac_colors[c("TLR7c_24", "BCRc_24", "DN2c_24")] <- plot_colors[c("TLR7c_24", "BCRc_48", "DN2c_48")] 

cluster_colors <- 
    c("C0" = "#A8CDE2", "C1" = "#3B83B9", "C2" = "#E3362C", "C3" = "#F9B56F", 
      "C4" = "#FC9230", "C5" = "#DDA086", "C6" = "#9F7BB8", "C7" = "#987898", 
      "C8" = "#F1E78D", "C9" = "#B05D2F", "C10" = "#83BF98", "C11" = "#6ABD5D", 
      "C12" = "#6F8544", "C13" = "#F4817F")
    

## Timecourse data
gene_list_bulk <- readr::read_lines("./data/bulk_genes.txt")

gene_exp <- readr::read_rds("./data/deseq_normalized_counts.rds")

## Gene tracks for ATAC-seq plot
AnnotationHub::setAnnotationHubOption("CACHE", "./data/AnnotationHub")
ah <- AnnotationHub::AnnotationHub(localHub=TRUE)
ens_data <- ah[["AH98047"]]

## Atac-seq tracks
bigwigs <- 
    "./data/bigwig" |>
    list.files(full.name = TRUE, pattern = "*.bigWig")

names(bigwigs) <- sub("^([^_]+_\\d+).+$", "\\1", basename(bigwigs))

## Cite-seq data
sc_cells <- 
    "./data/singlecell/cells.txt" |>
    readr::read_lines()

sc_genes <-
    "./data/singlecell/genes.txt" |>
    readr::read_lines()

sc_meta <- 
    "./data/singlecell/metadata.tsv" |>
    readr::read_tsv() |>
    dplyr::select(barcode, hto = dmm_hto_call, cluster = seurat_clusters) |>
    dplyr::mutate(cluster = factor(cluster, levels = paste0("C", 0:13)))

umap_df <-
    "./data/singlecell/umap_data.tsv" |>
    readr::read_tsv() |>
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

cluster_labels <-
    umap_df |>
    dplyr::group_by(cluster) |>
    dplyr::summarise_at(vars(UMAP_1, UMAP_2), mean) |>
    dplyr::ungroup()

sc_clusters_plot <- 
    ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(fill = cluster), 
               size = 1.5, 
               shape = 21, 
               stroke = .05, 
               color = "black") +
    geom_label(data = cluster_labels, 
               aes(x = UMAP_1, y = UMAP_2, label = cluster),
               label.padding = unit(0.1, "lines"),
               size = 12, size.unit = "pt", alpha = .5, fontface = "bold") +
    scale_fill_manual(values = cluster_colors) +
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
                
sc_hto_plot <- 
    ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(fill = hto), 
               size = 1.5, 
               shape = 21, 
               stroke = .05, 
               color = "black") +
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

sc_var_plots <- list("HTO" = sc_hto_plot, "Clusters" = sc_clusters_plot)

## Splicing
contrasts <- 
    c("Unstim 0h vs. TLR7c 24h" = "unstday0vs.TLR7",
      "Unstim 0h vs. BCRc 24h" = "unstday0vs.BCR",
      "Unstim 0h vs. DN2c 24h" = "unstday0vs.DN2",
      "TLR7c 24h vs. BCRc 24h" = "TLR7vs.BCR",
      "TLR7c 24h vs. DN2c 24h" = "TLR7vs.DN2",
      "BCRc 24h vs. DN2c 24h" = "BCRvs.DN2")


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
    navbarMenu("Single-cell RNA-seq", 
               tabPanel("UMAPs",
                        fluidRow(
                            column(6, selectizeInput("varsc", "Select variable:", choices = c("HTO", "Clusters"))),
                            column(6, selectizeInput("genesc", "Select gene:", choices = NULL))
                        ),
                        fluidRow(
                            column(6, plotOutput("plotscvars", width = "100%", height = "600px")),
                            column(6, plotOutput("plotsc", width = "100%", height = "600px"))
                        )
               ),
               tabPanel("Bubbleplot",
                        selectizeInput("markergenes", "Select genes:", multiple = TRUE, choices = NULL),
                        plotOutput("bubbleplot", inline = TRUE)
               )
    ),
    tabPanel("Splicing",
             selectizeInput("contrast", "Select contrast:", choices = names(contrasts)),
             fluidRow(
                 column(6,
                        div(id = "clusterTable",
                            h4(id = "title","Differential splicing events (clusters)"),
                            hr(),
                            div(DT::DTOutput("all_clusters"))
                        )
                 ),
                 column(6,
                        div(id="clusterView",
                            h4(id="title","Splicing event visualization"),
                            hr(),
                            div(plotOutput("select_cluster_plot", width = "100%")),
                            DT::DTOutput("cluster_view") 
                        )
                 )
             )
    )
)

server <- function(input, output, session) {

    ## Bulk RNA
    updateSelectizeInput(session, 
                         "generna", 
                         selected = "CD19",
                         choices = gene_list_bulk, 
                         server = TRUE)
    
    datarna <- reactive({
        req(input$generna)
        gene_exp |> format_timecourse(input$generna)
        })
    
    output$plotrna <- 
        renderPlot(res = 96, {
            req(input$generna)
            ggplot(data = datarna()) +
                       geom_quasirandom(aes(x = timep, y = norm_counts, fill = condition),
                                        size = 2.5, stroke = .25, shape = 21,
                                        method = "smiley", width = .2) +
                       scale_fill_manual(values = plot_colors) +
                       facet_grid(cols = vars(stim), space = "free", scale = "free_x") +
                       theme_bw() +
                       guides(fill = "none") +
                       labs(x = "hours", y = "Norm. counts")
            })
    
    ## Bulk ATAC
    updateSelectizeInput(session, 
                         "geneatac", 
                         selected = "CD19",
                         choices = gene_list_bulk, 
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
    output$plotscvars <- renderPlot(sc_var_plots[[input$varsc]])
    
    updateSelectizeInput(session, 
                         "genesc", 
                         selected = "CD19",
                         choices = sc_genes, 
                         server = TRUE)
    
    sc_data_gene <- 
        reactive({
            req(input$genesc)
            get_sc_data(input$genesc) |>
                tibble::add_column(barcode = sc_cells, .before = 1) |>
                dplyr::select(barcode, value = 2) |>
                dplyr::right_join(umap_df, dplyr::join_by(barcode)) |>
                dplyr::arrange(value) |> 
                dplyr::mutate(barcode = forcats::fct_inorder(barcode))
        })
    
    output$plotsc <- 
        renderPlot(
            ggplot(data = sc_data_gene(), 
                   aes(x = UMAP_1, y = UMAP_2)) +
                geom_point(aes(color = value), size = 1.5, stroke = 0) +
                scale_color_gradientn(colors = c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84",
                                                 "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000")) +
                theme_minimal() +
                theme(axis.text = element_text(size = 12),
                      axis.title = element_text(size = 12),
                      legend.text = element_text(size = 12),
                      legend.title = element_text(size = 12),
                      panel.grid = element_blank()) +
                guides(color = guide_colorbar(barwidth = .5, barheight = 8))
            )
    
    ### Bubble plot
    updateSelectizeInput(session, 
                         "markergenes", 
                         choices = sc_genes, 
                         server = TRUE)
    
    markers_plot_data <-
        reactive({
            req(input$markergenes)
            
            if (length(input$markergenes) == 1) {
                dat <- 
                    get_sc_data(input$markergenes) |>
                        tibble::add_column(barcode = sc_cells, .before = 1) |>
                        dplyr::left_join(sc_meta, dplyr::join_by(barcode)) |>
                        dplyr::select(-barcode, -hto) |>
                        tidyr::pivot_longer(-cluster, names_to = "gene") |>
                        dplyr::group_by(cluster, gene) |>
                        dplyr::summarise(avg.exp = mean(expm1(value)),
                                         pct.exp = mean(value > 0) * 100) |>
                        dplyr::ungroup()
            } else {
                dat <- 
                    purrr::map_dfc(input$markergenes, get_sc_data) |>
                        tibble::add_column(barcode = sc_cells, .before = 1) |>
                        dplyr::left_join(sc_meta, dplyr::join_by(barcode)) |>
                        dplyr::select(-barcode, -hto) |>
                        dplyr::select(cluster, dplyr::everything()) |>
                        dplyr::group_by(cluster) |>
                        dplyr::summarise_all(list(avg.exp = ~mean(expm1(.)),
                                                  pct.exp = ~mean(. > 0) * 100)) |>
                        dplyr::ungroup() |>
                        tidyr::pivot_longer(-cluster, 
                                            names_to = c("gene", ".value"),
                                            names_pattern = c("(.+)_(avg.exp|pct.exp)"))
                
                clust <- dat |>
                    dplyr::select(cluster, gene, avg.exp) |>
                    tidyr::pivot_wider(names_from = cluster, values_from = avg.exp) |>
                    tibble::column_to_rownames("gene") |>
                    dist() |>
                    hclust()
                
                dat <- dat |>
                    dplyr::mutate(gene = factor(gene, levels = clust$labels[clust$order]))
            }
            
            dat |>
                dplyr::group_by(gene) |>
                dplyr::mutate(avg.exp.scaled = as.numeric(scale(log1p(avg.exp))),
                              avg.exp.scaled = dplyr::case_when(avg.exp.scaled > 2.5 ~ 2.5,
                                                                avg.exp.scaled < -2.5 ~ -2.5,
                                                                .default = avg.exp.scaled)) |>
                dplyr::ungroup()
        })
    
    plot_height = reactive({
        req(input$markergenes)
        min_len <- 300
        custom_len <- 50 * length(input$markergenes)
        
        return(max(min_len, custom_len))
        })
    
    output$bubbleplot <- 
        renderPlot(res = 96, width = 800, height = function() plot_height(),
                   {
                       req(input$markergenes)
                       ggplot(markers_plot_data(), aes(x = cluster, y = gene)) +
                           geom_point(aes(size = pct.exp, fill = avg.exp.scaled),
                                      stroke = 0.2, shape = 21) +
                           scale_radius(range = c(0, 6),
                                        breaks = c(0, .25, .5, .75, 1) * 100,
                                        limits = c(0, 100)) +
                           scale_fill_gradient2(low = "Light Sky Blue", 
                                                mid = "lightyellow", 
                                                high = "Dark Red",
                                                midpoint = 0) +
                           theme_minimal() +
                           theme(
                               axis.text.x = element_text(size = 12),
                               axis.text.y = element_text(size = 12, face = 'italic'),
                               panel.grid = element_line(linewidth = .25, color = "grey90"),
                               legend.title = element_text(size = 12),
                               legend.key.spacing.y = unit(-.5, "lines")) +
                           guides(fill = guide_colorbar(order = 1, position = "right",
                                                        barwidth = .5, barheight = 5)) +
                           labs(x = NULL, y = NULL, fill = "Scaled\nExpression", size = "%\nExpressed")
                       }
                   )
    
    ## Splicing
    get_splicing_data <- 
        reactive({
            req(input$contrast)
            splicing_contrast <- contrasts[[input$contrast]]
            load(paste0("./data/splicing/", splicing_contrast, ".Rdata"))
            
            return(list("clusters" = clusters, 
                        "exons_table" = exons_table, 
                        "meta" = meta,
                        "cluster_ids" = cluster_ids,
                        "counts" = counts,
                        "introns" = introns))
        })
    
    splicing_data <- reactive({function() get_splicing_data()}())
    
    output$all_clusters <- 
        DT::renderDT({
            DT::datatable(splicing_data()$clusters[, c("gene", "coord", "N", "FDR", "annotation")],
                          escape = FALSE,
                          rownames = FALSE,
                          colnames = c("Genomic location" = "coord", "Gene" = "gene", "N" = "N", "Annotation" = "annotation", "q" = "FDR"),
                          selection = "single",
                          caption = "Click on a row to plot the corresponding visualization. N: number of introns within a cluster. q: Benjamini–Hochberg q-value.",
                          fillContainer = FALSE,
                          options = list(pageLength = 15,
                                         columnDefs = list(list(className = 'dt-center', targets = 0:4)))
            )
        })
    
    values <- reactiveValues(default = 1) 
    
    observeEvent(input$contrast, {values$default <- 1})
    
    observeEvent(input$all_clusters_rows_selected, {
        values$default <- input$all_clusters_rows_selected
    })
    
    plot_cluster_data <- 
        eventReactive(values$default, {
            sel <- values$default 
            gene  <- splicing_data()$clusters[ sel, ]$gene
            gene <- gsub("<.*?>", "", gene) # strip out html italic tags
            width <- leafviz::getGeneLength(splicing_data()$exons_table, gene)
            clusterID <- splicing_data()$clusters[ sel, ]$clusterID
            coord <- splicing_data()$clusters[ sel, ]$coord
            return(list(gene = gene, width = width, cluster = clusterID, coord = coord) )
    })
    
    output$select_cluster_plot <- 
        renderPlot(width = "auto", height = "auto", res = 90, {
            plotTitle <- c(plot_cluster_data()$gene, as.character(plot_cluster_data()$cluster) )
            leafviz::make_cluster_plot(plot_cluster_data()$cluster,
                                       main_title = plotTitle,
                                       meta = splicing_data()$meta,
                                       cluster_ids = splicing_data()$cluster_ids,
                                       exons_table = splicing_data()$exons_table,
                                       counts = splicing_data()$counts,
                                       introns = splicing_data()$introns)
        })
    
    output$cluster_view = DT::renderDT({
        clu <- plot_cluster_data()$cluster
        if(!is.null(clu)){
            if(length(splicing_data()$introns)){
                DT::datatable(leafviz::filter_intron_table(splicing_data()$introns, clu, toSave=FALSE),
                              rownames = TRUE,
                              options <- list(searching = FALSE, paging = FALSE, info = FALSE)
                )
            }
        } else {
            print("no cluster selected!")
        }
        
    }) 
    
    
}

shinyApp(ui, server)