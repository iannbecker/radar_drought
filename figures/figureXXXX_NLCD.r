##############################
#
#   Figure XXXX - NLCD response curves
#   Ian Becker
#   November 2025
#
##############################

library(mgcv)
library(gratia)  
library(ggplot2)
library(dplyr)

# Set working directo

# Configuration
season <- "fall"  # Change to "spring" as needed
output_dir <- paste0("gam_results/", season, "_model_comparison/nlcd_plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat(paste(rep("=", 60), collapse=""), "\n")
cat("GAM NLCD HABITAT PLOT -", toupper(season), "SEASON\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

#################################
# NLCD CODE LOOKUP TABLE
#################################

# NLCD classification with alpha codes
nlcd_lookup <- data.frame(
  code = c(11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 51, 52, 71, 81, 82, 90, 95),
  alpha_code = c("OW", "PI", "OSD", "LID", "MID", "HID", "BR", "DF", "EF", "MF", 
                 "DS", "SS", "GL", "PH", "CC", "WWo", "WWe"),
  name = c("Open Water", "Perennial Ice/Snow", "Developed, Open Space", 
           "Developed, Low Intensity", "Developed, Medium Intensity", 
           "Developed, High Intensity", "Barren Land", "Deciduous Forest", 
           "Evergreen Forest", "Mixed Forest", "Dwarf Scrub", "Shrub/Scrub", 
           "Grassland/Herbaceous", "Pasture/Hay", "Cultivated Crops", 
           "Woody Wetlands", "Emergent Herbaceous Wetlands"),
  stringsAsFactors = FALSE
)

cat("NLCD Lookup Table:\n")
print(nlcd_lookup)
cat("\n")

#################################
# LOAD MODEL
#################################

cat("Loading best model...\n")
model_file <- paste0("gam_results/", season, "_model_comparison/", season, "_best_model.rds")

if(!file.exists(model_file)) {
  stop("Model file not found: ", model_file)
}

best_model <- readRDS(model_file)
cat("Model loaded successfully\n\n")

#################################
# EXTRACT NLCD TERM
#################################

cat("Extracting NLCD smooth term...\n")
model_summary <- summary(best_model)

# Get all smooth terms
smooth_terms <- model_summary$s.table
smooth_names <- rownames(smooth_terms)

# Find NLCD interaction term
nlcd_term <- smooth_names[grepl("NLCD", smooth_names) & grepl("SPEI", smooth_names)]

if(length(nlcd_term) == 0) {
  stop("No NLCD interaction term found in model")
}

cat("Found NLCD term:", nlcd_term, "\n")
cat("F-statistic:", smooth_terms[nlcd_term, "F"], "\n")
cat("p-value:", smooth_terms[nlcd_term, "p-value"], "\n\n")

#################################
# GENERATE NLCD PLOT WITH ALPHA CODES
#################################

cat("Generating NLCD plot with alpha codes...\n")

tryCatch({
  # Extract smooth estimates using gratia
  nlcd_smooth <- gratia::smooth_estimates(best_model, smooth = nlcd_term, n = 100)
  
  cat("Smooth estimates extracted\n")
  cat("Columns in smooth estimates:", paste(names(nlcd_smooth), collapse=", "), "\n")
  
  # Identify the NLCD factor column
  nlcd_col <- names(nlcd_smooth)[grepl("NLCD", names(nlcd_smooth))][1]
  
  if(is.null(nlcd_col)) {
    stop("Could not identify NLCD column in smooth estimates")
  }
  
  cat("NLCD column identified:", nlcd_col, "\n")
  
  # Convert NLCD codes to alpha codes
  nlcd_smooth$nlcd_numeric <- as.numeric(as.character(nlcd_smooth[[nlcd_col]]))
  
  nlcd_smooth <- nlcd_smooth %>%
    left_join(nlcd_lookup, by = c("nlcd_numeric" = "code"))
  
  # Check for missing alpha codes
  if(any(is.na(nlcd_smooth$alpha_code))) {
    cat("WARNING: Some NLCD codes not found in lookup table\n")
    cat("Missing codes:", unique(nlcd_smooth$nlcd_numeric[is.na(nlcd_smooth$alpha_code)]), "\n")
  }
  
  # Create factor with alpha codes for plotting
  nlcd_smooth$nlcd_alpha <- factor(nlcd_smooth$alpha_code, 
                                   levels = nlcd_lookup$alpha_code)
  
  # Identify SPEI column
  spei_col <- names(nlcd_smooth)[grepl("SPEI", names(nlcd_smooth)) & 
                                   !grepl("by", names(nlcd_smooth))][1]
  
  cat("SPEI column identified:", spei_col, "\n")
  
  # Create the plot
  p <- ggplot(nlcd_smooth, aes(x = .data[[spei_col]], y = .estimate, group = nlcd_alpha)) +
    geom_line(aes(color = nlcd_alpha), linewidth = 1.2) +
    geom_ribbon(aes(ymin = .estimate - 2*.se, 
                    ymax = .estimate + 2*.se, 
                    fill = nlcd_alpha), 
                alpha = 0.2) +
    scale_color_viridis_d(option = "turbo", name = "Habitat Type") +
    scale_fill_viridis_d(option = "turbo", name = "Habitat Type") +
    theme_minimal(base_size = 14) +
    labs(title = paste(toupper(season), "Season: Drought Effects by Habitat Type"),
         subtitle = paste("Model term:", nlcd_term),
         x = "Drought Index (SPEI)",
         y = "Effect on Stopover Density") +
    theme(legend.position = "right",
          panel.grid.minor = element_blank())
  
  # Save plot
  plot_file <- file.path(output_dir, paste0(season, "_nlcd_drought_interaction_alpha.png"))
  ggsave(plot_file, p, width = 12, height = 8, dpi = 300)
  cat("NLCD plot saved:", plot_file, "\n")
  
  # Also create a faceted version for clarity
  p_facet <- ggplot(nlcd_smooth, aes(x = .data[[spei_col]], y = .estimate)) +
    geom_line(aes(color = nlcd_alpha), linewidth = 1) +
    geom_ribbon(aes(ymin = .estimate - 2*.se, 
                    ymax = .estimate + 2*.se, 
                    fill = nlcd_alpha), 
                alpha = 0.2) +
    facet_wrap(~nlcd_alpha, scales = "free_y", ncol = 4) +
    scale_color_viridis_d(option = "turbo", guide = "none") +
    scale_fill_viridis_d(option = "turbo", guide = "none") +
    theme_minimal(base_size = 12) +
    labs(title = paste(toupper(season), "Season: Drought Effects by Habitat Type (Faceted)"),
         subtitle = paste("Model term:", nlcd_term),
         x = "Drought Index (SPEI)",
         y = "Effect on Stopover Density") +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))
  
  facet_file <- file.path(output_dir, paste0(season, "_nlcd_drought_faceted_alpha.png"))
  ggsave(facet_file, p_facet, width = 14, height = 10, dpi = 300)
  cat("Faceted NLCD plot saved:", facet_file, "\n")
  
  # Save the data with alpha codes
  nlcd_data_file <- file.path(output_dir, paste0(season, "_nlcd_smooth_estimates_alpha.csv"))
  write.csv(nlcd_smooth, nlcd_data_file, row.names = FALSE)
  cat("NLCD smooth estimates saved:", nlcd_data_file, "\n")
  
}, error = function(e) {
  cat("ERROR generating NLCD plot:", e$message, "\n")
  traceback()
})

#################################
# SUMMARY
#################################

cat("\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("NLCD PLOT GENERATION COMPLETE\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("Output directory:", output_dir, "\n")
cat("Files created:\n")
cat("- Combined NLCD plot with alpha codes\n")
cat("- Faceted NLCD plot with alpha codes\n")
cat("- NLCD smooth estimates CSV with alpha codes\n")
cat(paste(rep("=", 60), collapse=""), "\n")