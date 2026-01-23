##############################
#
#   Figure 5 Panel C: Effect Size Comparison
#   Memory-efficient version (one model at a time)
#
##############################

library(mgcv)
library(gratia)
library(ggplot2)
library(ggpattern)
library(dplyr)

setwd("/Users/ianbecker/Library/CloudStorage/OneDrive-TheUniversityofTexas-RioGrandeValley/DroughtRadar/GAMs")

# Colors
fixed_colors <- c(
  "Concurrent Season" = "#7b3294",
  "Immediate Previous Season" = "#e08214"
)

# Categorization function
categorize_drought <- function(term) {
  if(grepl("1995.*2020|1995.*2024|1995.*2025|spring_1995", term)) {
    return("Concurrent Season")
  } else if(grepl("prev.*summer|prev_summer|prev.*winter|prev_winter", term)) {
    return("Immediate Previous Season")
  } else if(grepl("prev.*spring.*szn|prev_spring_szn|prev.*fall.*szn|prev_fall_szn", term)) {
    return("Previous Migratory Season")
  } else if(grepl("prev.*yr|prev.*year|prev_fall_yr|prev_spring_yr", term)) {
    return("Previous Year (1-year lag)")
  } else {
    return(NA)
  }
}

# Calculate effect size function
calculate_effect_size <- function(model, term, model_summary) {
  tryCatch({
    smooth_est <- smooth_estimates(model, select = term, n = 200)
    estimates <- smooth_est$.estimate
    se <- smooth_est$.se
    
    effect_range <- max(estimates) - min(estimates)
    idx_min <- which.min(estimates)
    idx_max <- which.max(estimates)
    se_range <- sqrt(se[idx_min]^2 + se[idx_max]^2)
    
    lower <- effect_range - 1.96 * se_range
    upper <- effect_range + 1.96 * se_range
    
    return(list(
      effect_size = effect_range,
      lower = lower,
      upper = upper,
      success = TRUE
    ))
  }, error = function(e) {
    return(list(success = FALSE, error = e$message))
  })
}

#################################
# PROCESS FALL MODEL
#################################

cat("Processing FALL model...\n")

# Load fall model
fall_model <- readRDS("fall_best_model.rds")
fall_summary <- summary(fall_model)

# Get drought terms
all_terms <- rownames(fall_summary$s.table)
drought_terms <- all_terms[grepl("SPEI", all_terms) & !grepl("ti\\(|te\\(|,", all_terms)]

# Create data frame
fall_data <- data.frame(
  term = drought_terms,
  f_stat = fall_summary$s.table[drought_terms, "F"],
  season = "Fall",
  stringsAsFactors = FALSE
)

fall_data$category <- sapply(fall_data$term, categorize_drought)
fall_data <- fall_data %>% filter(!is.na(category))

# Calculate effect sizes
for(i in 1:nrow(fall_data)) {
  result <- calculate_effect_size(fall_model, fall_data$term[i], fall_summary)
  if(result$success) {
    fall_data$effect_size[i] <- result$effect_size
    fall_data$lower[i] <- result$lower
    fall_data$upper[i] <- result$upper
  }
}

cat("  ✓ Fall effect sizes calculated\n")

# Save fall results and clear memory
saveRDS(fall_data, "fall_effect_sizes_temp.rds")
rm(fall_model, fall_summary)
gc()

cat("  ✓ Fall model cleared from memory\n\n")

#################################
# PROCESS SPRING MODEL
#################################

cat("Processing SPRING model...\n")

# Load spring model
spring_model <- readRDS("spring_best_model.rds")
spring_summary <- summary(spring_model)

# Get drought terms
all_terms <- rownames(spring_summary$s.table)
drought_terms <- all_terms[grepl("SPEI", all_terms) & !grepl("ti\\(|te\\(|,", all_terms)]

# Create data frame
spring_data <- data.frame(
  term = drought_terms,
  f_stat = spring_summary$s.table[drought_terms, "F"],
  season = "Spring",
  stringsAsFactors = FALSE
)

spring_data$category <- sapply(spring_data$term, categorize_drought)
spring_data <- spring_data %>% filter(!is.na(category))

# Calculate effect sizes
for(i in 1:nrow(spring_data)) {
  result <- calculate_effect_size(spring_model, spring_data$term[i], spring_summary)
  if(result$success) {
    spring_data$effect_size[i] <- result$effect_size
    spring_data$lower[i] <- result$lower
    spring_data$upper[i] <- result$upper
  }
}

cat("  ✓ Spring effect sizes calculated\n")

# Save spring results and clear memory
saveRDS(spring_data, "spring_effect_sizes_temp.rds")
rm(spring_model, spring_summary)
gc()

cat("  ✓ Spring model cleared from memory\n\n")

#################################
# COMBINE AND PLOT
#################################

cat("Loading saved results and creating plot...\n")

# Load both saved results
fall_data <- readRDS("fall_effect_sizes_temp.rds")
spring_data <- readRDS("spring_effect_sizes_temp.rds")

# Combine data
combined_data <- bind_rows(fall_data, spring_data)

combined_data <- combined_data %>%
  filter(category %in% c("Concurrent Season", "Immediate Previous Season"))

combined_data$category <- factor(
  combined_data$category,
  levels = c("Concurrent Season", "Immediate Previous Season")
)
combined_data$season <- factor(combined_data$season, levels = c("Fall", "Spring"))

# Create x-axis position - CATEGORY first, then season
combined_data$x_position <- factor(
  paste(combined_data$category, combined_data$season),
  levels = c(
    "Concurrent Season Fall",
    "Concurrent Season Spring",
    "Immediate Previous Season Fall",
    "Immediate Previous Season Spring"
  )
)

# Assign numeric positions - EVENLY SPACED (1, 2, 3, 4)
combined_data <- combined_data %>%
  mutate(x_numeric = case_when(
    x_position == "Concurrent Season Fall" ~ 1,
    x_position == "Concurrent Season Spring" ~ 2,
    x_position == "Immediate Previous Season Fall" ~ 3,
    x_position == "Immediate Previous Season Spring" ~ 4
  ))

# Create plot with legend
p <- ggplot(combined_data, aes(x = x_numeric, y = effect_size)) +
  geom_col_pattern(
    aes(fill = category, pattern = season),
    width = 0.75,
    color = "white",
    linewidth = 0.4,
    pattern_density = 0.3,
    pattern_spacing = 0.02,
    pattern_angle = 45,
    pattern_fill = "white",
    pattern_color = "gray30"
  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    linewidth = 0.8,
    color = "gray20"
  ) +
  scale_fill_manual(
    name = "Drought Timing",
    values = c("Concurrent Season" = "#7b3294", "Immediate Previous Season" = "#e08214")
  ) +
  scale_pattern_manual(
    values = c("Fall" = "none", "Spring" = "stripe"),
    guide = "none"  # Hide the pattern legend
  ) +
  scale_x_continuous(
    breaks = c(1, 2, 3, 4),
    labels = NULL,
    limits = c(0.4, 4.6),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = NULL,
    y = "Effect Size on Stopover"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 11),
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.margin = margin(5, 8, 5, 8),
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.9, "lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 40, 10)
  ) +
  guides(
    fill = guide_legend(override.aes = list(pattern = "none"))  # Force solid fills in legend
  )

ggsave("panel_c_effect_size.png", p, width = 8, height = 5, dpi = 300)

cat("  ✓ Plot saved\n\n")

# Clean up temp files
file.remove("fall_effect_sizes_temp.rds")
file.remove("spring_effect_sizes_temp.rds")

cat("Complete!\n")


