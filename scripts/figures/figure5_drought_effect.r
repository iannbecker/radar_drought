##############################
#
#   Figure 5 - Individual Drought Effects
#   Ian Becker
#   December 2025
#
#   Generates individual plots for each drought smooth term
#
##############################

library(mgcv)
library(gratia)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/Users/ianbecker/Library/CloudStorage/OneDrive-TheUniversityofTexas-RioGrandeValley/DroughtRadar/GAMs")

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
  
  # Create individual plot for each term
  for(term in drought_terms) {
    
    cat("  Plotting:", term, "...\n")
    
    tryCatch({
      # Extract smooth estimates
      smooth_est <- smooth_estimates(current_model, select = term, n = 200)
      
      # Get SPEI column name
      spei_col <- names(smooth_est)[grepl("SPEI", names(smooth_est))][1]
      
      # Extract data
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
      
      # Create plot
      p <- ggplot(plot_data, aes(x = SPEI, y = estimate)) +
        geom_hline(yintercept = 0, linetype = "dashed", 
                   color = "gray40", linewidth = 0.5) +
        geom_ribbon(aes(ymin = lower, ymax = upper), 
                    fill = plot_color, alpha = 0.3) +
        geom_line(color = plot_color, linewidth = 1.2) +
        labs(
          title = paste(toupper(season), "-", category),
          subtitle = paste("F =", round(f_stat, 0)),
          x = "SPEI (Drought Index)",
          y = "Effect on Stopover (log scale)"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 12, color = "gray30"),
          panel.grid.minor = element_blank(),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(size = 12)
        )
      
      # Save plot
      clean_term <- gsub("[^[:alnum:]]", "_", term)
      filename <- file.path(output_dir, paste0(season, "_", clean_term, ".png"))
      ggsave(filename, p, width = 8, height = 6, dpi = 300)
      
      cat("    ✓ Saved:", basename(filename), "\n")
      
    }, error = function(e) {
      cat("    ✗ Error:", e$message, "\n")
    })
  }
  
  cat("\n")
}

cat("===============================================\n")
cat("ALL INDIVIDUAL PLOTS COMPLETE!\n")
cat("===============================================\n")
cat("Output directory:", output_dir, "\n")
cat("\nFiles organized by:\n")
cat("  - fall_*.png: Fall drought smooths\n")
cat("  - spring_*.png: Spring drought smooths\n")
cat("\nEach plot shows:\n")
cat("  - Season and drought timing category in title\n")
cat("  - F-statistic in subtitle\n")
cat("  - Color-coded by drought timing\n")
cat("  - 95% confidence intervals (ribbons)\n")