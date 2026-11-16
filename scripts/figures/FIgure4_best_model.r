##############################
#
#   Figure 4 - Effect Size Plots
#   Ian Becker
#   November 2026
#
##############################

# Plots each smooth term's effect size (range of its partial effect on the
# log scale, i.e. peak-to-trough) from the {season}_best_model_terms.csv
# tables written by results_pipeline.R. Interaction ranges there are already
# restricted to data-supported bins.

library(ggplot2)
library(dplyr)

# Settings ------------------------------------

results_dir <- "/Users/ianbecker/Desktop/project_code/radar_drought/data/final_analysis/results_outputs"      # folder holding {season}_best_model_terms.csv
season      <- "fall"         # "fall" or "spring"
plot_height <- 14               # inches; increase for more space between labels

# Terms included as controls, not plotted
exclude_terms <- c("s(year)", "s(dist_radar_800m)")

# DATA PREP  ------------------------------------

terms <- read.csv(file.path(results_dir, paste0(season, "_best_model_terms.csv")))

terms <- terms %>%
  filter(!term %in% exclude_terms, !is.na(effect_range))

# Readable labels for every variable name in either season's model

var_labels <- c(
  SPEI_fall_current         = "SPEI current fall",
  SPEI_spring_prev_year     = "SPEI previous spring",
  temp_fall_PRISM_stack     = "Temperature",
  temp_spring_PRISM_stack   = "Temperature",
  ppt_fall_PRISM_stack      = "Precipitation",
  ppt_spring_PRISM_stack    = "Precipitation",
  fall_uwind                = "U-wind",
  uwind_spring              = "U-wind",
  fall_vwind                = "V-wind",
  vwind_spring              = "V-wind",
  msavi_fall                = "Vegetation greenness",
  msavi_spring              = "Vegetation greenness",
  alan_dist_stack           = "Distance to artificial light",
  dwater_perm_stack         = "Distance to permanent water",
  dwater_temp_stack         = "Distance to temporary water",
  dcoast                    = "Distance to coast",
  landscape_pct_wetland     = "% Wetland",
  landscape_pct_agriculture = "% Agriculture",
  landscape_pct_developed   = "% Developed",
  landscape_pct_grassland   = "% Grassland",
  landscape_pct_shrub_scrub = "% Shrub/scrub",
  landscape_pct_forest      = "% Forest"
)

climate_vars <- c("temp_fall_PRISM_stack", "temp_spring_PRISM_stack",
                  "ppt_fall_PRISM_stack", "ppt_spring_PRISM_stack",
                  "fall_uwind", "uwind_spring", "fall_vwind", "vwind_spring")

# Split "s(x)" / "ti(x,y)" into variable names
term_vars <- function(term) {
  strsplit(gsub("^(s|ti)\\(|\\)$", "", term), ",")[[1]]
}

terms$vars <- lapply(terms$term, term_vars)

terms$term_clean <- sapply(terms$vars, function(v) {
  lab <- var_labels[v]
  if (any(is.na(lab))) stop("No label for: ", paste(v[is.na(lab)], collapse = ", "))
  paste(lab, collapse = " × ")
})


terms$category <- sapply(terms$vars, function(v) {
  if (length(v) > 1)          return("Interaction")
  if (grepl("^SPEI", v))      return("Drought")
  if (v %in% climate_vars)    return("Climate")
  "Habitat"
})

# Order largest effect at top

terms <- terms %>% arrange(desc(effect_range))
terms$term_clean <- factor(terms$term_clean, levels = rev(terms$term_clean))

# PLOTTING ------------------------------------

category_colors <- c(
  "Drought"     = "#c94c4c",
  "Climate"     = "#4c7ac9",
  "Habitat"     = "#52a354",
  "Interaction" = "#9370DB"
)

p <- ggplot(terms, aes(x = effect_range, y = term_clean, color = category)) +
  geom_segment(aes(x = 0, xend = effect_range, yend = term_clean),
               linewidth = 0.7, alpha = 0.6) +
  geom_point(size = 3.5, alpha = 0.9) +
  scale_color_manual(values = category_colors, name = "Term Type") +
  scale_x_continuous(
    name = "Effect size (range of partial effect, log scale)",
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10)
  )

print(p)

# Save plot (label spacing is set by the height)

ggsave(file.path(results_dir, paste0(season, "_effect_size_plot.png")),
       plot = p, width = 10, height = plot_height, dpi = 300, bg = "white")
cat("\n✓ Effect size plot saved!\n")