library(tidyverse)
library(ggbeeswarm)

fig_colors <-
    "./figure_colors.txt" |>
    read_tsv(col_names = c("stim", "timep", "col")) |>
    filter(timep == 24,
	   stim %in% c("IL-4c", "TLR7c", "BCRc", "DN2c")) |>
    select(stim, col) |>
    deframe()

ldsc_df <- 
    "../atacseq/ldsc/compiled_results.tsv" |>
    read_tsv() |>
    mutate(set = recode(set, 
			"IL4" = "IL-4c", 
			"TLR7" = "TLR7c", 
			"BCR" = "BCRc",
			"DN2" = "DN2c"),
	   set = factor(set, levels = c("IL-4c", "TLR7c", "BCRc", "DN2c")),
	   group = case_when(grepl("T2D|Height|HEIGHTz|Type_2_Diabetes|Schizophrenia|MDD|LDL|HDL|cancer", trait) ~ "control",
			     TRUE ~ "test")) |>
    filter(group == "test") |>
    mutate(group = case_when(grepl("Lupus|Multiple_scle|Rheumatoid_Arthritis|Crohns|_CD_", trait) ~ "Group 1",
			     grepl("Type_1_Dia|IBD|Celiac|_UC_|Ulcerative_Col|Primary_bil|PSORIASIS", trait) ~ "Group 2",
			     grepl("Asthma|ASTHMA|ALLERGY", trait) ~ "Group 3",
			     .default = "Other"),
	   group = factor(group, levels = paste("Group", 1:3)))

ldsc_summ_df <-
    ldsc_df |>
    mutate(trait2 = case_when(grepl("ASTHMA|GBMI\\.Asthma", trait) ~ "asthma",
                              grepl("AdultOnsetAsthma", trait) ~ "adult_asthma",
                              grepl("ChildOnsetAsthma", trait) ~ "child_asthma",
                              grepl("Celiac", trait) ~ "celiac",
                              grepl("Covid", trait) ~ "covid",
                              grepl("Crohns|_CD_", trait) ~ "crohns",
                              grepl("ALLERGY", trait) ~ "allergy",
                              grepl("_IBD", trait) ~ "ibd",
                              grepl("Multiple_sclerosis", trait) ~ "ms",
                              grepl("biliary", trait) ~ "pbc",
                              grepl("PSORIASIS", trait) ~ "psoriasis",
                              grepl("Rheumatoid", trait) ~ "ra",
                              grepl("Lupus", trait) ~ "sle",
                              grepl("Type_1_Diabetes", trait) ~ "t1d",
                              grepl("Ulcerative|_UC_", trait) ~ "uc")) |>
    group_by(set, trait2) |>
    summarize(tau_star = mean(tau_star))

plot_out <-
    ggplot(ldsc_summ_df, aes(x = set, y = tau_star)) +
    geom_boxplot(aes(fill = set), outlier.shape = NA, alpha = .25) +
	geom_quasirandom(aes(fill = set), 
			 method = "smiley", 
			 shape = 21, stroke = .2, size = 1.5, width = .3) +
	#scale_fill_manual(values = c("#2C3E50", "#D35400", "#95A5A6")) +
	scale_fill_manual(values = fig_colors) +
    theme_minimal() +
	theme(panel.grid.minor.y = element_blank(),
	      axis.text = element_text(size = 10),
	      axis.title = element_text(size = 12),
	      plot.background = element_rect(fil = "white", color = "white")) +
	labs(x = NULL, y = expression(tau * "*")) +
    guides(fill = "none")

ggsave("./ldsc_tau.png", plot_out, width = 2.75, height = 3)
