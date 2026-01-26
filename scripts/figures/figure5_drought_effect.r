##############################
#
#   Figure 5 - Individual Drought Effects
#   Ian Becker
#   December 2025
#
#
##############################

# This script is used to create all of the components for Figure 5
# in the main body of the paper showing effect size of drought between
# seasons and individual drought effects 

library(mgcv)
library(gratia)
library(ggplot2)
library(ggpattern)
library(dplyr)
library(tidyr)

setwd("PATH TO MODELS HERE")

#################################
# PREP
#################################

# Output directory

output_dir <- "individual_drought_plots"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load in models 

cat("Loading models...\n")
fall_model <- readRDS("fall_best_model.rds")
spring_model <- readRDS("spring_best_model.rds")

fall_summary <- summary(fall_model)
spring_summary <- summary(spring_model)
cat("Models loaded successfully\n\n")

# Create function to categorize drought terms

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

# Color scheme

fixed_colors <- c(
  "Concurrent Season" = "#7b3294",
  "Immediate Previous Season" = "#e08214",
  "Previous Migratory Season" = "#1b9e77",
  "Previous Year (1-year lag)" = "#377eb8"
)

#################################
# INDIVIDUAL DROUGHT PLOTS
#################################

cat("CREATING INDIVIDUAL DROUGHT SMOOTH PLOTS...\n\n")

for(season in c("spring", "fall")) {
  
  cat("Processing", toupper(season), "season...\n")
  
  # Select appropriate model and summary
  
  current_model <- if(season == "fall") fall_model else spring_model
  current_summary <- if(season == "fall") fall_summary else spring_summary
  
  # Get drought terms
  
  all_terms <- rownames(current_summary$s.table)
  drought_terms <- all_terms[grepl("SPEI", all_terms) & 
                               !grepl("ti\\(|te\\(|,", all_terms)]
  
  cat("Found", length(drought_terms), "drought smooth terms\n")
  
  # Create individual plot for each drought window
  
  for(term in drought_terms) {
    
    cat("  Plotting:", term, "...\n")
    
    tryCatch({
      
      # Extract smooth estimates
      
      smooth_est <- smooth_estimates(current_model, select = term, n = 200)
      
      # Get SPEI column
      
      spei_col <- names(smooth_est)[grepl("SPEI", names(smooth_est))][1]
      
      # Create plot data
      
      plot_data <- data.frame(
        SPEI = smooth_est[[spei_col]],
        estimate = smooth_est$.estimate,
        se = smooth_est$.se,
        lower = smooth_est$.estimate - 1.96 * smooth_est$.se,
        upper = smooth_est$.estimate + 1.96 * smooth_est$.se
      )
      
      # Get category and color
      
      category <- categorize_drought(term)
      if(is.na(category)) category <- "Other"
      
      plot_color <- if(category %in% names(fixed_colors)) {
        fixed_colors[[category]]
      } else {
        "gray50"
      }
      
      # Get F-statistic
      
      f_stat <- current_summary$s.table[term, "F"]
      
      # Get axis ranges for breaks
      
      x_range <- range(plot_data$SPEI)
      y_range <- range(c(plot_data$lower, plot_data$upper))
      
      # Create plot with minimal axes
      
      p <- ggplot(plot_data, aes(x = SPEI, y = estimate)) +
        geom_hline(yintercept = 0, linetype = "dashed", 
                   color = "gray40", linewidth = 0.3) +
        geom_ribbon(aes(ymin = lower, ymax = upper), 
                    fill = plot_color, alpha = 0.3) +
        geom_line(color = plot_color, linewidth = 1.5) +
        scale_x_continuous(
          breaks = c(x_range[1], 0, x_range[2]),
          labels = c(round(x_range[1], 1), "0", round(x_range[2], 1))
        ) +
        scale_y_continuous(
          breaks = c(y_range[1], 0, y_range[2]),
          labels = c(round(y_range[1], 1), "0", round(y_range[2], 1))
        ) +
        labs(x = NULL, y = NULL) +
        theme_minimal(base_size = 20) +  
        theme(
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
          axis.text = element_text(size = 18, face = "bold"),  
          axis.ticks = element_line(linewidth = 0.5),
          axis.ticks.length = unit(0.15, "cm"),
          plot.margin = margin(5, 5, 5, 5)
        )
      
      # Save plot
      
      clean_term <- gsub("[^[:alnum:]]", "_", term)
      filename <- file.path(output_dir, paste0(season, "_", clean_term, "_small.png"))
      ggsave(filename, p, width = 4, height = 3, dpi = 300)
      
      cat("    ✓ Saved:", basename(filename), "\n")
      
    }, error = function(e) {
      cat("    ✗ Error:", e$message, "\n")
    })
  }
}

cat("===============================================\n")
cat("ALL INDIVIDUAL PLOTS COMPLETE!\n")
cat("===============================================\n")

#################################
# BAR CHART - EFFECT SIZES
#################################

cat("CREATING EFFECT SIZE BAR CHART...\n\n")

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

# Process FALL

cat("Processing FALL effect sizes...\n")

all_terms <- rownames(fall_summary$s.table)
drought_terms <- all_terms[grepl("SPEI", all_terms) & !grepl("ti\\(|te\\(|,", all_terms)]

fall_data <- data.frame(
  term = drought_terms,
  f_stat = fall_summary$s.table[drought_terms, "F"],
  season = "Fall",
  stringsAsFactors = FALSE
)

fall_data$category <- sapply(fall_data$term, categorize_drought)
fall_data <- fall_data %>% filter(!is.na(category))

for(i in 1:nrow(fall_data)) {
  result <- calculate_effect_size(fall_model, fall_data$term[i], fall_summary)
  if(result$success) {
    fall_data$effect_size[i] <- result$effect_size
    fall_data$lower[i] <- result$lower
    fall_data$upper[i] <- result$upper
  }
}

cat("  ✓ Fall effect sizes calculated\n")

# Clear fall model

rm(fall_model, fall_summary)
gc()

# Process SPRING

cat("Processing SPRING effect sizes...\n")

all_terms <- rownames(spring_summary$s.table)
drought_terms <- all_terms[grepl("SPEI", all_terms) & !grepl("ti\\(|te\\(|,", all_terms)]

spring_data <- data.frame(
  term = drought_terms,
  f_stat = spring_summary$s.table[drought_terms, "F"],
  season = "Spring",
  stringsAsFactors = FALSE
)

spring_data$category <- sapply(spring_data$term, categorize_drought)
spring_data <- spring_data %>% filter(!is.na(category))

for(i in 1:nrow(spring_data)) {
  result <- calculate_effect_size(spring_model, spring_data$term[i], spring_summary)
  if(result$success) {
    spring_data$effect_size[i] <- result$effect_size
    spring_data$lower[i] <- result$lower
    spring_data$upper[i] <- result$upper
  }
}

cat("  ✓ Spring effect sizes calculated\n")

# Clear spring model

rm(spring_model, spring_summary)
gc()

# Combine data

combined_data <- bind_rows(fall_data, spring_data)

combined_data <- combined_data %>%
  filter(category %in% c("Concurrent Season", "Immediate Previous Season"))

combined_data$category <- factor(
  combined_data$category,
  levels = c("Concurrent Season", "Immediate Previous Season")
)
combined_data$season <- factor(combined_data$season, levels = c("Fall", "Spring"))

# Create x-axis position

combined_data$x_position <- factor(
  paste(combined_data$category, combined_data$season),
  levels = c(
    "Concurrent Season Fall",
    "Concurrent Season Spring",
    "Immediate Previous Season Fall",
    "Immediate Previous Season Spring"
  )
)

# Assign numeric positions

combined_data <- combined_data %>%
  mutate(x_numeric = case_when(
    x_position == "Concurrent Season Fall" ~ 1,
    x_position == "Concurrent Season Spring" ~ 2,
    x_position == "Immediate Previous Season Fall" ~ 3,
    x_position == "Immediate Previous Season Spring" ~ 4
  ))

# Create bar chart

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
    guide = "none"
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
    fill = guide_legend(override.aes = list(pattern = "none"))
  )

ggsave(file.path(output_dir, "effect_size_bar_chart.png"), p, width = 8, height = 5, dpi = 300)

cat("===============================================\n")
cat("BAR CHART SAVED!\n")
cat("===============================================\n")
