##############################
#
#   Quick Smooth Curve Viewer
#   Ian Becker
#   December 2025
#
##############################

# Use this to quickly view any of the smooth curves from GAM output

library(mgcv)
library(gratia)
library(ggplot2)
library(dplyr)

setwd("/Users/ianbecker/Library/CloudStorage/OneDrive-TheUniversityofTexas-RioGrandeValley/DroughtRadar/GAMs")

#################################
# CONFIGURATION
#################################

# Choose season

season <- "spring"

# List the predictors you want to plot (exact term names from model)

predictors_to_plot <- c(
  "s(dwater_perm_stack)",
  "s(dwater_temp_stack)",
  "s(alan_dist_stack)"
)

# Create readable labels for predictors

predictor_labels <- c(
  "s(dwater_perm_stack)" = "Distance to Permanent Water",
  "s(dwater_temp_stack)" = "Distance to Temporary Water",
  "s(alan_dist_stack)" = "Distance to ALAN"
)

#################################
# LOAD MODEL
#################################

cat("Loading", season, "model...\n")
model_file <- paste0(season, "_best_model.rds")
model <- readRDS(model_file)
model_summary <- summary(model)

cat("Model loaded successfully\n\n")

#################################
# CHECK AVAILABLE TERMS
#################################

cat("Available smooth terms in model:\n")
all_terms <- rownames(model_summary$s.table)
for(term in all_terms) {
  f_stat <- model_summary$s.table[term, "F"]
  p_val <- model_summary$s.table[term, "p-value"]
  sig <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "")))
  cat("  ", term, " (F =", round(f_stat, 1), sig, ")\n")
}
cat("\n")

#################################
# EXTRACT SMOOTH EFFECTS
#################################

# Get the actual model data
model_data <- model$model

cat("=== PLOTTING SMOOTHS ===\n\n")

for(term in predictors_to_plot) {
  
  var_name <- gsub("s\\(|\\)", "", term)
  
  cat("Processing:", var_name, "\n")
  
  if(var_name %in% names(model_data)) {
    
    # Get actual data range
    actual_range <- range(model_data[[var_name]], na.rm = TRUE)
    actual_mean <- mean(model_data[[var_name]], na.rm = TRUE)
    actual_sd <- sd(model_data[[var_name]], na.rm = TRUE)
    
    cat("  Actual data range:", actual_range[1], "to", actual_range[2], "\n")
    cat("  Mean:", actual_mean, ", SD:", actual_sd, "\n")
    
    # METHOD 1: Let gratia handle everything (EASIEST)
    # Just specify which smooth, let it figure out the rest
    smooth_est <- smooth_estimates(
      model, 
      select = term,
      n = 200
      # Don't pass data argument - let it use model$model
    )
    
    pred_col <- names(smooth_est)[grepl(var_name, names(smooth_est))][1]
    
    if(is.null(pred_col)) {
      pred_col <- names(smooth_est)[2]  # Fallback to second column
    }
    
    # Get F-statistic
    f_stat <- model_summary$s.table[term, "F"]
    p_val <- model_summary$s.table[term, "p-value"]
    
    sig_text <- ifelse(p_val < 0.001, "p < 0.001", 
                       ifelse(p_val < 0.01, "p < 0.01",
                              ifelse(p_val < 0.05, "p < 0.05",
                                     paste("p =", round(p_val, 3)))))
    
    # Create plot
    p <- ggplot(smooth_est, aes(x = .data[[pred_col]], y = .estimate)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
      geom_ribbon(aes(ymin = .estimate - 1.96 * .se, 
                      ymax = .estimate + 1.96 * .se), 
                  alpha = 0.3, fill = "#377eb8") +
      geom_line(color = "#377eb8", linewidth = 1) +
      # Add rug plot to show data distribution
      geom_rug(data = model_data[sample(nrow(model_data), min(1000, nrow(model_data))), ],
               aes(x = .data[[var_name]], y = NULL), 
               alpha = 0.2, sides = "b", color = "#377eb8", length = unit(0.02, "npc")) +
      labs(
        title = var_name,
        subtitle = paste0("F = ", round(f_stat, 1), ", ", sig_text),
        x = paste0(var_name, " (standardized)"),
        y = "Effect on Stopover Density (log scale)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "gray30"),
        axis.title = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank()
      )
    
    filename <- paste0(season, "_", var_name, "_standardized.png")
    ggsave(filename, p, width = 8, height = 6, dpi = 300)
    cat("  Saved:", filename, "\n")
    
    # Also print the actual x-axis range to verify
    cat("  Plot x-axis range:", range(smooth_est[[pred_col]]), "\n")
    
  } else {
    cat("  WARNING: Variable", var_name, "not found in model data!\n")
  }
  
  cat("\n")
}

cat("DONE!\n")
