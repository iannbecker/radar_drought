##############################
#
#   Figure XXXX - Forest Plots
#   Ian Becker
#   November 2025
#
##############################

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

# Clean up term names for the supplementary table

clean_term_names <- function(term) {
  term <- gsub("s\\(", "", term)
  term <- gsub("ti\\(", "", term)
  term <- gsub("\\)", "", term)
  term <- gsub("_stack", "", term)
  term <- gsub("_fall_PRISM", "", term)
  term <- gsub("_PRISM", "", term)
  term <- gsub("SPEI_1995_2020", "SPEI current season", term)
  term <- gsub("SPEI_prev_fall_yr", "SPEI previous fall year", term)
  term <- gsub("SPEI_prev_spring_szn", "SPEI previous spring", term)
  term <- gsub("SPEI_prev_summer", "SPEI previous summer", term)
  term <- gsub(",", " × ", term)
  term <- gsub("alan_dist", "Distance to artificial light", term)
  term <- gsub("dwater_perm", "Distance to permanent water", term)
  term <- gsub("dwater_temp", "Distance to temporary water", term)
  term <- gsub("forest_age", "Forest age", term)
  term <- gsub("temp_fall", "Temperature", term)
  term <- gsub("ppt_fall", "Precipitation", term)
  term <- gsub("wetness_fall", "Soil wetness", term)
  term <- gsub("fall_uwind", "U-wind", term)
  term <- gsub("fall_vwind", "V-wind", term)
  term <- gsub("msavi_fall", "Vegetation greenness", term)
  term <- gsub("dcoast", "Distance to coast", term)
  return(term)
}

smooth_df$term_full <- clean_term_names(smooth_df$term)

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

smooth_df$category <- sapply(smooth_df$term_full, assign_category)

# Calculate F-statistics and confidence intervals

smooth_df$F_stat <- smooth_df$F

# Calculate approximate 95% CI 

alpha <- 0.05
smooth_df$F_lower <- smooth_df$F_stat * qchisq(alpha/2, df = smooth_df$Ref.df) / smooth_df$Ref.df
smooth_df$F_upper <- smooth_df$F_stat * qchisq(1 - alpha/2, df = smooth_df$Ref.df) / smooth_df$Ref.df

# Sort by F-statistic (largest to smallest)

smooth_df <- smooth_df %>%
  arrange(desc(F_stat)) %>%
  mutate(term_number = row_number())

# Create numbered labels for the plot

smooth_df$term_label <- paste0("Term ", smooth_df$term_number)

# Reorder factor levels for plotting (top to bottom)

smooth_df$term_label <- factor(smooth_df$term_label, 
                               levels = rev(smooth_df$term_label))

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

p <- ggplot(smooth_df, aes(x = F_stat, y = term_label, color = category)) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.5) + # Reference line at 0
  geom_errorbarh(aes(xmin = F_lower, xmax = F_upper),
                 height = 0.3, linewidth = 0.7, alpha = 0.7) + # Error bars
  geom_point(size = 3.5, alpha = 0.9) + # Points
  scale_color_manual(values = category_colors, name = "Term Type") + # Use color scale 
  scale_x_continuous(labels = scales::comma) + # Format x-axis labels
  labs(x = "F-statistic", y = NULL) +
  theme_bw(base_size = 15) +
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
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  )

print(p) 

# Save the plot

ggsave(file.path(model_dir, "spring_forest_plot.png"), plot = p, width = 10, height = 7, dpi = 300, bg = "white")
cat("\n✓ Forest plot saved!\n")

################
#  Supplementary Table
################

# Create supplementary table with full term names

supp_table <- smooth_df %>%
  arrange(term_number) %>%
  select(
    Term_Number = term_number,
    Term_Name = term_full,
    Category = category,
    F_Statistic = F_stat,
    EDF = edf,
    Ref_DF = Ref.df,
    p_value = `p-value`
  ) %>%
  mutate(
    F_Statistic = round(F_Statistic, 2),
    EDF = round(EDF, 2),
    Ref_DF = round(Ref_DF, 2),
    p_value = ifelse(p_value < 0.001, "< 0.001", format.pval(p_value, digits = 3))
  )

# Save supplementary table

write.xlsx(supp_table, 
          file.path(model_dir, "spring_forest_plot_supplementary_table.xlsx"),
          rowNames = FALSE)

cat("✓ Supplementary table saved!\n")

# Print table

cat("=== SUPPLEMENTARY TABLE ===\n")
print(supp_table)


cat("\n✓ Figure XXXX complete!\n")
