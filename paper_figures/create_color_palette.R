library(tidyverse)
library(colorspace)

# The optimal 16-point step size to preserve hues at 72h
lum_scale <- c(
  "0"  = 98,  # Background Gray
  "4"  = 88,  # Pastel/Light
  "24" = 72,  # Medium-Light
  "48" = 56,  # Medium-Dark
  "72" = 40   # Dark (Floor)
)

hues <- c(
  "Unstim"    = NA,  
  "BCR/TLR7c" = 10,  # Rust / Brick
  "CD40c"     = 50,  # Orange
  "IL-4c"     = 85,  # Yellow/Gold
  "TLR7c"     = 120, # Warm Leaf Green (Pulled far away from blue)
  "BCRc"      = 250, # Royal Blue / Indigo (Pulled far away from green)
  "TLR9c"     = 290, # Purple
  "DN2c"      = 330  # Pink/Magenta
)

palette_df <- data.frame(Condition = character(), Time = character(), Hex = character())

for (cond in names(hues)) {
  if (cond == "Unstim") {
    t_points <- c("0", "4", "24")
  } else if (cond == "IL-4c") {
    t_points <- c("4", "24", "72") 
  } else {
    t_points <- c("4", "24", "48", "72")
  }
  
  for (t in t_points) {
    L <- lum_scale[t]
    
    if (cond == "Unstim") {
      hex_code <- hex(polarLUV(L = L, C = 0, H = 0))
    } else {
      H <- hues[cond]
      # Enforce strict L matching
      max_c <- max_chroma(h = H, l = L)
      C <- min(80, max_c) 
      hex_code <- hex(polarLUV(L = L, C = C, H = H))
    }
    palette_df <- rbind(palette_df, data.frame(Condition = cond, Time = t, Hex = hex_code))
  }
}

palette_df <- 
    palette_df |>
    mutate(Condition = factor(Condition, levels = c("Unstim", "IL-4c", "CD40c", "TLR9c", "TLR7c", "BCRc", "BCR/TLR7c", "DN2c")), 
	   Time = factor(Time, levels = c(0, 4, 24, 48, 72))) |>
    unite("lab", c(Condition, Time), sep = "_", remove = FALSE) |>
    arrange(Condition, Time) |>
    as_tibble()

print(palette_df, n = Inf)

# Plot
p <- 
    ggplot(palette_df, 
	   aes(x = Condition, y = Time)) +
    geom_point(aes(color = lab), size = 12) +
    scale_color_manual(values = set_names(palette_df$Hex, palette_df$lab)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(color = "none")

ggsave("testcolors.png", p, width = 7, height = 7)

write_tsv(palette_df, "./figure_colors_final.txt")
