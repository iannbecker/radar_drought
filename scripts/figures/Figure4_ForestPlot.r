##############################
#
#   Figure 4 - Forest Plots
#   Ian Becker
#   November 2025
#
##############################

# This script makes figure 4 in the main text

library(mgcv)
library(ggplot2)
library(dplyr)
library(openxlsx)

# Set working directory and load model

model_dir <- "PATH HERE"
best_model <- readRDS(file.path(model_dir, "spring_best_model.rds"))

################
#  Data Prep
################

# Extract smooth term statistics

smooth_summary <- summary(best_model)$s.table
smooth_df <- as.data.frame(smooth_summary)
smooth_df$term <- rownames(smooth_df)

# Filter out terms we don't want

exclude_patterns <- c(
  ":ecoregion", 
  ":NLCD", 
  "^s\\(ecoregion\\)",
  "^s\\(year\\)"
)

keep_terms <- !grepl(paste(exclude_patterns, collapse = "|"), smooth_df$term)
smooth_df <- smooth_df[keep_terms, ]

# Clean up term names for the plot axes

clean_term_names <- function(term) {
  term <- gsub("s\\(", "", term)
  term <- gsub("ti\\(", "", term)
  term <- gsub("\\)", "", term)
  term <- gsub("_stack", "", term)
  term <- gsub("SPEI_spring_1995_2020", "SPEI current season", term)
  term <- gsub("SPEI_prev_fall_szn", "SPEI previous fall", term)
  term <- gsub("SPEI_prev_spring_yr", "SPEI previous spring year", term)
  term <- gsub("SPEI_prev_winter", "SPEI previous winter", term)
  term <- gsub(",", " × ", term)
  term <- gsub("alan_dist", "Distance to artificial light", term)
  term <- gsub("dwater_perm", "Distance to permanent water", term)
  term <- gsub("dwater_temp", "Distance to temporary water", term)
  term <- gsub("forest_age", "Forest age", term)
  term <- gsub("temp_spring_PRISM", "Temperature", term)
  term <- gsub("ppt_spring_PRISM", "Precipitation", term)
  term <- gsub("wetness_spring", "Soil wetness", term)
  term <- gsub("uwind_spring", "U-wind", term)
  term <- gsub("vwind_spring", "V-wind", term)
  term <- gsub("NLCD_1995_2020", "NLCD Land Cover", term)
  term <- gsub("msavi_spring", "Vegetation greenness", term)
  term <- gsub("radar_id", "Radar", term)
  term <- gsub("dcoast", "Distance to coast", term)
  return(term)
}

smooth_df$term_clean <- clean_term_names(smooth_df$term)

# Assign categories

assign_category <- function(term) {
  term_lower <- tolower(term)
  
  # Check if it's an interaction
  if (grepl(" × ", term)) {
    return("Interaction")
  }
  
  # Main effects
  if (grepl("spei", term_lower)) return("Drought")
  if (grepl("temp|ppt|wetness|wind|soil", term_lower)) return("Climate")
  if (grepl("light|coast|water|forest|vegetation|msavi", term_lower)) return("Habitat")
  
  return("Other")
}

smooth_df$category <- sapply(smooth_df$term_clean, assign_category)

# Calculate F-statistics and confidence intervals

smooth_df$F_stat <- smooth_df$F

# Calculate approximate 95% CI 
alpha <- 0.05
smooth_df$F_lower <- smooth_df$F_stat * qchisq(alpha/2, df = smooth_df$Ref.df) / smooth_df$Ref.df
smooth_df$F_upper <- smooth_df$F_stat * qchisq(1 - alpha/2, df = smooth_df$Ref.df) / smooth_df$Ref.df

# Sort by F-statistic (largest to smallest)
smooth_df <- smooth_df %>%
  arrange(desc(F_stat))

# Reorder factor levels for plotting (top to bottom, largest F first)
smooth_df$term_clean <- factor(smooth_df$term_clean, 
                               levels = rev(smooth_df$term_clean))

cat("\n✓ Data prep complete, time to plot!\n")

################
#  Plotting
################

# Define color scheme

category_colors <- c(
  "Drought" = "#c94c4c",
  "Climate" = "#4c7ac9", 
  "Habitat" = "#52a354",
  "Interaction" = "#9370DB"
)

# Create the forest plot

p <- ggplot(smooth_df, aes(x = F_stat, y = term_clean, color = category)) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = F_lower, xmax = F_upper),
                 height = 0.3, linewidth = 0.7, alpha = 0.7) +
  geom_point(size = 3.5, alpha = 0.9) +
  scale_color_manual(values = category_colors, name = "Term Type") +
  scale_x_continuous(labels = scales::comma) +
  labs(x = "F-statistic", y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray30"),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 0),
    axis.text.y = element_text(size = 9),  # Adjust size for readability
    axis.text.x = element_text(size = 10)
  )

print(p) 

# Save the plot (make it taller to accommodate term names)

ggsave(file.path(model_dir, "spring_forest_plot.png"), 
       plot = p, width = 10, height = 12, dpi = 300, bg = "white")
cat("\n✓ Forest plot saved!\n")

cat("\n✓ Figure 4 complete!\n")