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

model_dir <- "PATH HERE"

#################################
# Data prep
#################################

# NLCD alpha code table

nlcd_lookup <- data.frame(
  code = c(11, 12, 21, 22, 23, 24, 31, 41, 42, 43, 51, 52, 71, 81, 82, 90, 95),
  alpha_code = c("OWat", "PI", "DevOS", "DevLI", "DevMI", "DevHI", "BarL", "DecidF", "EverF", "MixF", 
                 "DwaSc", "SH/SC", "Gra/Herb", "Past/Hay", "CulCro", "WoodWE", "EmergWe"),
  name = c("Open Water", "Perennial Ice/Snow", "Developed, Open Space", 
           "Developed, Low Intensity", "Developed, Medium Intensity", 
           "Developed, High Intensity", "Barren Land", "Deciduous Forest", 
           "Evergreen Forest", "Mixed Forest", "Dwarf Scrub", "Shrub/Scrub", 
           "Grassland/Herbaceous", "Pasture/Hay", "Cultivated Crops", 
           "Woody Wetlands", "Emergent Herbaceous Wetlands"),
  stringsAsFactors = FALSE
)

cat("NLCD Alpha Codes:\n")
print(nlcd_lookup)

# Load model

best_model <- readRDS(file.path(model_dir, "spring_best_model.rds")) # change based on season

# Extract NLCD smooth terms

cat("Model Summary...\n")
model_summary <- summary(best_model)
print(model_summary)

# Get all smooth terms

smooth_terms <- model_summary$s.table
smooth_names <- rownames(smooth_terms)

# Find NLCD interaction term

nlcd_term <- smooth_names[grepl("NLCD", smooth_names) & grepl("SPEI", smooth_names)]

cat("Found NLCD term:", nlcd_term, "\n")
cat("F-statistic:", smooth_terms[nlcd_term, "F"], "\n")
cat("Extracting predictions from the model...\n")

# Get the model data

model_data <- best_model$model

# Identify SPEI and NLCD variables

spei_var <- names(model_data)[grepl("SPEI.*1995.*2020", names(model_data))][1]
nlcd_var <- names(model_data)[grepl("NLCD", names(model_data))][1]

cat("SPEI variable:", spei_var, "\n")
cat("NLCD variable:", nlcd_var, "\n")

# Get unique NLCD levels

nlcd_levels <- levels(model_data[[nlcd_var]])
cat("NLCD levels found:", paste(nlcd_levels, collapse=", "), "\n")

# Create prediction grid

spei_range <- seq(min(model_data[[spei_var]], na.rm = TRUE),
                  max(model_data[[spei_var]], na.rm = TRUE),
                  length.out = 100)

# Create prediction data for each NLCD level

pred_list <- list()

for(nlcd_level in nlcd_levels) {
  # Create prediction data frame matching model structure
  pred_data <- data.frame(
    spei = spei_range,
    nlcd = factor(nlcd_level, levels = nlcd_levels)
  )
  names(pred_data) <- c(spei_var, nlcd_var)
  
  # Add other variables at their median values
  for(var in names(model_data)) {
    if(!var %in% c(spei_var, nlcd_var, "stopover")) {
      if(is.numeric(model_data[[var]])) {
        pred_data[[var]] <- median(model_data[[var]], na.rm = TRUE)
      } else if(is.factor(model_data[[var]])) {
        pred_data[[var]] <- factor(levels(model_data[[var]])[1], 
                                   levels = levels(model_data[[var]]))
      }
    }
  }
  
  # Get predictions with standard errors
  preds <- predict(best_model, newdata = pred_data, se.fit = TRUE, 
                   type = "terms", terms = nlcd_term)
  
  # Store results
  pred_list[[nlcd_level]] <- data.frame(
    spei = spei_range,
    nlcd_code = nlcd_level,
    estimate = preds$fit[, 1],
    se = preds$se.fit[, 1]
  )
}

# Combine all predictions

nlcd_smooth <- do.call(rbind, pred_list)
names(nlcd_smooth)[1] <- spei_var

cat("Predictions extracted successfully\n")
cat("Columns in smooth estimates:", paste(names(nlcd_smooth), collapse=", "), "\n")

# Identify the NLCD factor column

nlcd_col <- "nlcd_code"
spei_col <- spei_var

cat("NLCD column identified:", nlcd_col, "\n")
cat("SPEI column identified:", spei_col, "\n")

# Convert NLCD codes to alpha codes

nlcd_smooth$nlcd_numeric <- as.numeric(as.character(nlcd_smooth[[nlcd_col]]))

nlcd_smooth <- nlcd_smooth %>%
  left_join(nlcd_lookup, by = c("nlcd_numeric" = "code"))

# Create factor with alpha codes for plotting

nlcd_smooth$nlcd_alpha <- factor(nlcd_smooth$alpha_code, 
                                 levels = nlcd_lookup$alpha_code)

# Remove Perennial Ice/Snow if it still exists

nlcd_smooth <- nlcd_smooth %>%
  filter(alpha_code != "PI")

#################################
#  Generate NLCD plots
################################

# Create the plot

p <- ggplot(nlcd_smooth, aes(x = .data[[spei_col]], y = estimate, group = nlcd_alpha)) +
  geom_line(aes(color = nlcd_alpha), linewidth = 1.2) +
  scale_color_hue(name = "Land Cover Type") +
  theme_minimal(base_size = 14) +
  labs(x = "SPEI", y = "Effect on Stopover Density") +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

# Save plot

plot_file <- file.path(model_dir, "spring_nlcd_drought_interaction_alpha.png")
ggsave(plot_file, p, width = 12, height = 8, dpi = 300)
cat("NLCD plot saved:", plot_file, "\n")

# Faceted version for supplementary

p_facet <- ggplot(nlcd_smooth, aes(x = .data[[spei_col]], y = estimate)) +
  geom_line(aes(color = nlcd_alpha), linewidth = 1) +
  geom_ribbon(aes(ymin = estimate - 2*se, 
                  ymax = estimate + 2*se, 
                  fill = nlcd_alpha), 
              alpha = 0.2) +
  facet_wrap(~nlcd_alpha, scales = "free_y", ncol = 4) +
  scale_color_hue(name = "Habitat Type") +
  theme_minimal(base_size = 12) +
  labs(x = "SPEI", y = "Effect on Stopover Density") +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

# Save supplementary faceted plot

facet_file <- file.path(model_dir, "spring_nlcd_drought_faceted_alpha.png")
ggsave(facet_file, p_facet, width = 14, height = 10, dpi = 300)
cat("Faceted NLCD plot saved:", facet_file, "\n")

# Save the data with alpha codes

nlcd_data_file <- file.path(model_dir, "spring_nlcd_smooth_estimates_alpha.csv")
write.csv(nlcd_smooth, nlcd_data_file, row.names = FALSE)
cat("NLCD smooth estimates saved:", nlcd_data_file, "\n")
