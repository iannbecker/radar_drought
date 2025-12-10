##############################
#
#   Figure 5 - Drought effects
#   Ian Becker
#   December 2025
#
##############################

library(mgcv)
library(gratia)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("PATH HERE")

#################################
# PREP
#################################

# Load in models 

cat("Loading models...\n")
fall_model <- readRDS("gam_results/fall_model_comparison/fall_best_model.rds")
spring_model <- readRDS("gam_results/spring_model_comparison/spring_best_model.rds")

# Check that models loaded in correctly

fall_summary <- summary(fall_model)
spring_summary <- summary(spring_model)
cat("Models loaded successfully\n\n")

# Create function to categorize drought terms (will be used later)

categorize_drought <- function(term) {
  if(grepl("1995.*2020|1995.*2024|1995.*2025", term)) {
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

# Color scheme (for later)

fixed_colors <- c(
  "Concurrent Season" = "#7b3294",
  "Immediate Previous Season" = "#e08214",
  "Previous Migratory Season" = "#1b9e77",
  "Previous Year (1-year lag)" = "#377eb8"
)

#################################
# Fall and Spring drought effects
#################################

cat("PLOTTING FALL AND SPRING DROUGHT EFFECTS...\n\n")

for(season in c("fall", "spring")) {
  
  cat("Processing", toupper(season), "season...\n")
  
  # Select appropriate model and summary
  
  current_model <- if(season == "fall") fall_model else spring_model
  current_summary <- if(season == "fall") fall_summary else spring_summary
  
  cat("Identifying drought smooth terms...\n")
  all_terms <- rownames(current_summary$s.table)
  
  drought_terms <- all_terms[grepl("SPEI", all_terms) & 
                               !grepl("ti\\(|te\\(|,", all_terms)]
  
  # Confirming all drought terms are accounted for
  
  cat("Found", length(drought_terms), "drought smooth terms:\n")
  for(term in drought_terms) {
    f_stat <- current_summary$s.table[term, "F"]
    cat("  -", term, "(F =", round(f_stat, 0), ")\n")
  }
  
  cat("\nExtracting drought effects...\n")
  
  all_drought_data <- list()
  
  # Processing drought smooth effects
  
  for(term in drought_terms) {
    cat("  Processing:", term, "...\n")
    
    tryCatch({
      smooth_est <- smooth_estimates(current_model, select = term, n = 200)
      
      spei_col <- names(smooth_est)[grepl("SPEI", names(smooth_est))][1]
      x_vals <- smooth_est[[spei_col]]
      
      y_vals <- smooth_est$.estimate
      se_vals <- smooth_est$.se
      
      ymin_vals <- y_vals - 1.96 * se_vals
      ymax_vals <- y_vals + 1.96 * se_vals
      
      f_stat <- current_summary$s.table[term, "F"]
      
      clean_label <- gsub("s\\(|\\)", "", term)
      readable_label <- categorize_drought(term)
      if(is.na(readable_label)) readable_label <- clean_label
      
      all_drought_data[[term]] <- data.frame(
        SPEI = x_vals,
        estimate = y_vals,
        se = se_vals,
        lower = ymin_vals,
        upper = ymax_vals,
        term = term,
        label = readable_label,
        f_stat = f_stat,
        stringsAsFactors = FALSE
      )
      
      cat("    Success - extracted", length(x_vals), "points\n")
      
    }, error = function(e) {
      cat("    ERROR processing", term, ":", e$message, "\n")
    })
  }
  
  drought_combined <- bind_rows(all_drought_data)
  
  drought_combined <- drought_combined %>%
    mutate(
      label_ordered = reorder(label, -f_stat)
    )
  
  cat("\nDrought data extracted!\n")
  
  cat("\nCreating plot...\n")
  
  drought_labels <- levels(drought_combined$label_ordered)
  color_mapping <- sapply(drought_labels, function(lbl) {
    if(as.character(lbl) %in% names(fixed_colors)) {
      fixed_colors[[as.character(lbl)]]
    } else {
      "gray50"
    }
  })
  names(color_mapping) <- drought_labels
  
  ribbon_alpha <- 0.4
  
  p1 <- ggplot(drought_combined, aes(x = SPEI, y = estimate, 
                                     color = label_ordered, 
                                     fill = label_ordered)) +
    geom_hline(yintercept = 0, linetype = "dashed", 
               color = "gray40", linewidth = 0.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                alpha = ribbon_alpha, color = NA) +
    geom_line(linewidth = .4, alpha = 0.9) +
    scale_color_manual(
      values = color_mapping,
      name = "Drought Timing",
      breaks = c("Concurrent Season", 
                 "Immediate Previous Season",
                 "Previous Migratory Season", 
                 "Previous Year (1-year lag)")
    ) +
    scale_fill_manual(
      values = color_mapping,
      name = "Drought Timing",
      breaks = c("Concurrent Season", 
                 "Immediate Previous Season",
                 "Previous Migratory Season", 
                 "Previous Year (1-year lag)")
    ) +
    labs(
      x = "SPEI (Standardized Precipitation-Evapotranspiration Index)",
      y = "Effect on Stopover Density (log scale)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray30"),
      plot.caption = element_text(size = 9, color = "gray50", hjust = 0),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10)
    )
  
  output_file_1 <- paste0(season, "_all_drought_smooths.png")
  ggsave(output_file_1, p1, width = 12, height = 8, dpi = 300)
  
  cat("\n", toupper(season), "plot saved:", output_file_1, "\n\n")
}

#################################
# Fall vs. Spring effect size
#################################

cat("\nDROUGHT EFFECT SIZE COMPARISON\n\n")

cat("Calculating effect sizes...\n")

# Function to calculate drought effect size for both seasons

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
    
    mean_abs_effect <- mean(abs(estimates))
    
    return(list(
      effect_size = effect_range,
      lower = lower,
      upper = upper,
      mean_abs = mean_abs_effect,
      success = TRUE
    ))
  }, error = function(e) {
    return(list(success = FALSE, error = e$message))
  })
}

# Extract effect sizes through both seasons

all_season_data <- list()

for(season in c("Fall", "Spring")) {
  
  cat("Processing", season, "effect sizes...\n")
  
  # Select appropriate model and summary
  
  current_model <- if(season == "Fall") fall_model else spring_model
  current_summary <- if(season == "Fall") fall_summary else spring_summary
  
  # Get drought terms
  
  all_terms <- rownames(current_summary$s.table)
  drought_terms <- all_terms[grepl("SPEI", all_terms) & !grepl("ti\\(|te\\(|,", all_terms)]
  
  # Create data frame
  
  season_data <- data.frame(
    term = drought_terms,
    f_stat = current_summary$s.table[drought_terms, "F"],
    season = season,
    stringsAsFactors = FALSE
  )
  
  # Categorize
  
  season_data$category <- sapply(season_data$term, categorize_drought)
  season_data <- season_data %>% filter(!is.na(category))
  
  # Calculate effect sizes
  
  for(i in 1:nrow(season_data)) {
    result <- calculate_effect_size(current_model, season_data$term[i], current_summary)
    if(result$success) {
      season_data$effect_size[i] <- result$effect_size
      season_data$lower[i] <- result$lower
      season_data$upper[i] <- result$upper
      season_data$mean_abs[i] <- result$mean_abs
    }
  }
  
  all_season_data[[season]] <- season_data
  cat(season, "effect sizes calculated\n")
}

cat("\n")

# Combine both seasons

combined_data <- bind_rows(all_season_data)

combined_data$category <- factor(
  combined_data$category,
  levels = c("Concurrent Season", 
             "Immediate Previous Season",
             "Previous Migratory Season",
             "Previous Year (1-year lag)")
)

combined_data$season <- factor(combined_data$season, levels = c("Fall", "Spring"))

combined_data <- combined_data %>%
  filter(category %in% c("Concurrent Season", "Immediate Previous Season"))

combined_data$category <- factor(
  combined_data$category,
  levels = c("Concurrent Season", "Immediate Previous Season")
)

cat("Effect sizes extracted:\n")
print(combined_data %>% select(category, season, effect_size, f_stat))

# Plot!

p2 <- ggplot(combined_data, aes(x = category, y = effect_size, fill = category, alpha = season)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "white",
    linewidth = 0.4
  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.3,
    linewidth = 0.8,
    color = "gray20"
  ) +
  scale_fill_manual(values = fixed_colors, guide = "none") +
  scale_alpha_manual(
    values = c("Fall" = 1, "Spring" = 0.6),
    name = "Season",
    guide = guide_legend(override.aes = list(fill = "gray50"))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = NULL,
    y = "Effect on Stopover (log scale)",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, color = "gray30", hjust = 0.5),
    plot.caption = element_text(size = 8, color = "gray50", hjust = 0),
    axis.text.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

output_file_2 <- "panel_c_effect_size.png"
ggsave(output_file_2, p2, width = 8, height = 5, dpi = 300)

cat("\nPart 2 plot saved:", output_file_2, "\n")

cat("\nALL DONE!\n")
