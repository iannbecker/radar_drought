##############################
#
#   GAM Response Curves - Manual Selection
#   Ian Becker
#   12/22/2025
#
##############################

library(mgcv)
library(gratia)
library(ggplot2)

setwd("/Users/ianbecker/Library/CloudStorage/OneDrive-TheUniversityofTexas-RioGrandeValley/DroughtRadar/GAMs")

#################################
# CONFIGURATION
#################################

season <- "spring"  # Change to "spring" as needed

# >>> SPECIFY WHICH TERMS TO PLOT <<<
terms_to_plot <- c("s(dwater_perm_stack)", "s(dwater_temp_stack)", "s(alan_dist_stack)")

# >>> SCALING OPTIONS <<<
# Fixed y-axis limits for all univariate plots (NULL = auto)
fixed_y_limits <- NULL  # Set to c(-0.5, 0.5) for fixed range

# Number of breaks on x-axis (for consistency)
n_x_breaks <- 5

# Plot dimensions (width x height in inches)
plot_width <- 8
plot_height <- 6

# Paths

model_dir <- getwd()
output_dir <- model_dir

#################################
# LOAD MODEL AND PREPROCESSING
#################################

cat("Loading model and preprocessing info...\n")

# Load model
model_file <- file.path(model_dir, paste0(season, "_best_model.rds"))
best_model <- readRDS(model_file)
cat("✓ Model loaded\n")

# Load preprocessing info
preprocessing_file <- file.path(model_dir, paste0(season, "_preprocessing_info.rds"))
preprocessing_info <- readRDS(preprocessing_file)
cat("✓ Preprocessing info loaded\n")

# Extract original stats
original_stats <- preprocessing_info$original_stats
numeric_vars <- preprocessing_info$numeric_vars
spatial_vars <- preprocessing_info$spatial_vars

cat("  Numeric vars (standardized):", length(numeric_vars), "\n")
cat("  Spatial vars (centered):", length(spatial_vars), "\n\n")

#################################
# HELPER FUNCTION: BACK-TRANSFORM
#################################

back_transform <- function(var_name, standardized_values) {
  if(var_name %in% names(original_stats)) {
    stats <- original_stats[[var_name]]
    
    if(var_name %in% numeric_vars) {
      # Standardized: back-transform using mean and SD
      return(standardized_values * stats$sd + stats$mean)
    } else if(var_name %in% spatial_vars) {
      # Centered: back-transform using mean only
      return(standardized_values + stats$mean)
    }
  }
  # No transformation needed
  return(standardized_values)
}

#################################
# PLOT EACH TERM
#################################

for(term in terms_to_plot) {
  
  cat("Plotting:", term, "\n")
  
  # Extract variable name
  var_name <- gsub("s\\(|ti\\(|te\\(|\\)|,.*", "", term)
  var_name <- trimws(var_name)
  
  # Get smooth estimates
  sm <- smooth_estimates(best_model, select = term, n = 200)
  
  # Check if this is a univariate smooth that needs back-transformation
  if(var_name %in% names(original_stats) && !grepl(",", term)) {
    
    cat("  Back-transforming", var_name, "to original scale...\n")
    
    # Get x-axis range in standardized units
    x_range <- range(sm[[var_name]], na.rm = TRUE)
    standardized_breaks <- seq(x_range[1], x_range[2], length.out = n_x_breaks)
    
    # Back-transform breaks
    original_breaks <- back_transform(var_name, standardized_breaks)
    
    # Create nice round breaks based on range
    range_size <- max(original_breaks) - min(original_breaks)
    
    if(range_size > 10000) {
      # Round to nearest 5000
      clean_labels <- round(original_breaks / 5000) * 5000
    } else if(range_size > 1000) {
      # Round to nearest 1000
      clean_labels <- round(original_breaks / 1000) * 1000
    } else if(range_size > 100) {
      # Round to nearest 100
      clean_labels <- round(original_breaks / 100) * 100
    } else if(range_size > 10) {
      # Round to nearest 10
      clean_labels <- round(original_breaks / 10) * 10
    } else if(range_size > 1) {
      # Round to nearest 1
      clean_labels <- round(original_breaks)
    } else {
      # Round to 1 decimal for small ranges
      clean_labels <- round(original_breaks, 1)
    }
    
    # Plot with back-transformed axis
    p <- ggplot(sm, aes(x = .data[[var_name]], y = .estimate)) +
      geom_ribbon(aes(ymin = .estimate - 2*.se, ymax = .estimate + 2*.se),
                  fill = "gray80", alpha = 0.5) +
      geom_line(color = "#2166ac", linewidth = 1.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      scale_x_continuous(
        name = var_name,
        breaks = standardized_breaks,
        labels = clean_labels
      ) +
      labs(
        y = "Partial effect on stopover",
        title = term
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12)
      )
    
    # Apply fixed y-limits if specified
    if(!is.null(fixed_y_limits)) {
      p <- p + coord_cartesian(ylim = fixed_y_limits)
    }
    
  } else {
    # Use gratia::draw for interactions or non-transformed variables
    p <- draw(best_model, select = term) +
      theme_minimal(base_size = 14) +
      theme(panel.grid.minor = element_blank())
  }
  
  # Save with consistent dimensions
  clean_name <- gsub("[^[:alnum:]]", "_", term)
  filename <- file.path(output_dir, paste0(season, "_", clean_name, ".png"))
  ggsave(filename, p, width = plot_width, height = plot_height, dpi = 300)
  
  cat("  ✓ Saved:", basename(filename), "\n\n")
}

cat("Done! Plots saved to:", output_dir, "\n")
