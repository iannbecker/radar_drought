##############################
#
#   GAMS Fall - OPTIMIZED WITH AUDIT
#   Ian Becker
#   9/25/2025 
#   
#
##############################

library(terra)
library(mgcv)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)

# Set working directory
setwd("/home/ianbecker01/drought_radar")

# Configuration
season <- "spring"  
sample_size <- 50000

# Input directories  
data_dir <- paste0(season, "_raster")  
output_dir <- paste0("gam_results/", season, "_model_comparison")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Memory monitoring function
monitor_memory <- function(stage) {
  cat("Memory at", stage, ":\n")
  print(gc())
  cat("---\n")
}

# BCR lookup table
bcr_lookup <- data.frame(
  code = c(0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
  name = c(
    "Central Hardwoods",
    "Central Mixed Grass Prairie",
    "Chihuahuan Desert",
    "Eastern Tallgrass Prairie",
    "Edwards Plateau",
    "Gulf Coastal Prairie",
    "Mississippi Alluvial Valley",
    "Oaks And Prairies",
    "Planicie Costera",
    "Shortgrass Prairie",
    "Sierra Madre Occidental",
    "Sierra Madre Oriental",
    "Southeastern Coastal Plain",
    "Southern Rockies/Colorado Plateau",
    "Tamaulipan Brushlands",
    "West Gulf Coastal Plain/Ouachitas"
  )
)

# Eco table printing function
print_eco_table <- function(eco_values, label = "") {
  cat(label, "\n")
  
  # Handle both numeric and factor inputs
  if(is.factor(eco_values)) {
    eco_codes <- as.numeric(levels(eco_values))[eco_values]
  } else {
    eco_codes <- as.numeric(eco_values)
  }
  
  eco_table <- table(eco_codes, useNA = "always")
  eco_df <- data.frame(
    BCR_code = as.numeric(names(eco_table)),
    count = as.numeric(eco_table)
  )
  eco_df <- merge(eco_df, bcr_lookup, by.x = "BCR_code", by.y = "code", all.x = TRUE)
  eco_df <- eco_df[order(eco_df$BCR_code), ]
  
  print(eco_df[, c("BCR_code", "name", "count")])
  cat("\n")
}

cat(paste(rep("=", 60), collapse=""), "\n")
cat("GAM MODEL COMPARISON -", toupper(season), "SEASON\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# Create audit log
audit_log <- list()

#################################
# DATA PREPARATION
#################################

cat("Loading raster stacks...\n")

# Load stopover data
stopover_file <- file.path(data_dir, paste0(season, "_stack.tif"))
if(!file.exists(stopover_file)) {
  alt_files <- list.files(data_dir, pattern = "stopover|radar", full.names = TRUE)
  if(length(alt_files) > 0) {
    stopover_file <- alt_files[1]
    cat("Using stopover file:", basename(stopover_file), "\n")
  } else {
    stop("Cannot find stopover data file in ", data_dir)
  }
}
stopover_stack <- rast(stopover_file)

cat("Stopover data: ", nlyr(stopover_stack), "layers (years)\n")
years <- 1995:(1995 + nlyr(stopover_stack) - 1)
cat("Analysis years:", min(years), "to", max(years), "\n")

# Load predictor variables
predictor_files <- list.files(data_dir, pattern = "\\.tif$", full.names = TRUE)
predictor_files <- predictor_files[!grepl("stopover|fall_stack|spring_stack", basename(predictor_files))]

predictor_stacks <- list()
for(file in predictor_files) {
  name <- tools::file_path_sans_ext(basename(file))
  name <- gsub(paste0("_", season, "_cropped"), "", name) 
  
  cat("Loading predictor:", name, "\n")
  raster_stack <- rast(file)
  predictor_stacks[[name]] <- raster_stack
  cat("  Loaded:", name, "with", nlyr(raster_stack), "layers\n")
}

cat("\nTotal predictor stacks:", length(predictor_stacks), "\n")

# AUDIT 1: Load ecoregion raster and check
cat("\n=== AUDIT 1: Ecoregion Raster ===\n")
ecoregion_raster <- predictor_stacks[["ecoregion"]]
audit_log$raster_unique_values <- sort(unique(values(ecoregion_raster)))
cat("Unique BCR codes in ecoregion raster:\n")
print(audit_log$raster_unique_values)

# Sample pixels for analysis
cat("\nSampling pixels for model training...\n")
set.seed(123)

# Use pixels with valid stopover in ANY year
any_valid_stopover <- app(stopover_stack, fun = function(x) any(!is.na(x)))
valid_pixels <- which(values(any_valid_stopover) == 1)
cat("Total valid pixels (non-NA stopover in any year):", length(valid_pixels), "\n")

# AUDIT 2: Check ecoregions in valid pixels
cat("\n=== AUDIT 2: Ecoregions in Valid Pixels ===\n")
audit_log$valid_pixel_ecoregions <- table(values(ecoregion_raster)[valid_pixels])
print_eco_table(values(ecoregion_raster)[valid_pixels], "BCR codes in valid pixels:")

# Sample pixels
sample_size <- min(sample_size, length(valid_pixels))
sampled_pixels <- sample(valid_pixels, sample_size)

cat("Sample size for training:", length(sampled_pixels), "\n")
cat("Efficiency: Using", round(length(sampled_pixels)/length(valid_pixels)*100, 1), "% of valid pixels\n")

# AUDIT 3: Check ecoregions in sampled pixels
cat("\n=== AUDIT 3: Ecoregions in Sampled Pixels ===\n")
audit_log$sampled_pixel_ecoregions <- table(values(ecoregion_raster)[sampled_pixels])
print_eco_table(values(ecoregion_raster)[sampled_pixels], "BCR codes in sampled pixels:")

# Create long-form dataset
cat("\nCreating long-form dataset...\n")
start_time <- Sys.time()

# Extract all predictor values
cat("Extracting all predictor values...\n")
predictor_data <- list()
n_years <- nlyr(stopover_stack)

for(pred_name in names(predictor_stacks)) {
  cat("  Extracting values for:", pred_name, "\n")
  pred_stack <- predictor_stacks[[pred_name]]
  
  if(nlyr(pred_stack) == 1) {
    static_vals <- values(pred_stack)[sampled_pixels]
    predictor_data[[pred_name]] <- rep(static_vals, n_years)
    cat("    Static predictor - replicated across", n_years, "years\n")
  } else if(nlyr(pred_stack) == n_years) {
    time_vals <- values(pred_stack)[sampled_pixels, ]
    predictor_data[[pred_name]] <- as.vector(t(time_vals))
    cat("    Time-varying predictor -", nlyr(pred_stack), "layers matched\n")
  } else {
    layer_names <- names(pred_stack)
    layer_years <- regmatches(layer_names, regexpr("\\d{4}", layer_names))
    
    if(length(layer_years) > 0) {
      matched_values <- matrix(NA, nrow = length(sampled_pixels), ncol = n_years)
      for(year_idx in 1:n_years) {
        current_year <- years[year_idx]
        matching_layer <- which(layer_years == current_year)
        if(length(matching_layer) > 0) {
          matched_values[, year_idx] <- values(pred_stack[[matching_layer[1]]])[sampled_pixels]
        }
      }
      predictor_data[[pred_name]] <- as.vector(t(matched_values))
      cat("    Multi-layer predictor - matched by year\n")
    } else {
      if(nlyr(pred_stack) >= n_years) {
        time_vals <- values(pred_stack[[1:n_years]])[sampled_pixels, ]
        predictor_data[[pred_name]] <- as.vector(t(time_vals))
        cat("    Multi-layer predictor - used first", n_years, "layers\n")
      } else {
        for(layer_idx in 1:nlyr(pred_stack)) {
          layer_name <- paste0(pred_name, "_L", layer_idx)
          static_vals <- values(pred_stack[[layer_idx]])[sampled_pixels]
          predictor_data[[layer_name]] <- rep(static_vals, n_years)
        }
        cat("    Created", nlyr(pred_stack), "separate predictors\n")
      }
    }
  }
  rm(pred_stack)
  gc()
}

# Extract stopover values
cat("Extracting stopover values...\n")
stopover_matrix <- values(stopover_stack)[sampled_pixels, ]
stopover_vector <- as.vector(t(stopover_matrix))

rm(stopover_stack, stopover_matrix)
gc()

# Create year and pixel ID vectors
pixel_ids <- rep(sampled_pixels, n_years)
year_vector <- rep(years, each = length(sampled_pixels))

# Build final dataset
cat("Assembling final dataset...\n")
pooled_data <- data.frame(
  pixel_id = pixel_ids,
  stopover = stopover_vector,
  year = year_vector,
  stringsAsFactors = FALSE
)

for(pred_name in names(predictor_data)) {
  pooled_data[[pred_name]] <- predictor_data[[pred_name]]
}

# AUDIT 4: Initial pooled_data
cat("\n=== AUDIT 4: Initial pooled_data ===\n")
cat("Unique ecoregion values:", sort(unique(pooled_data$ecoregion)), "\n")
audit_log$initial_pooled <- table(pooled_data$ecoregion)
print_eco_table(pooled_data$ecoregion, "BCR distribution in pooled_data:")
cat("Total observations:", nrow(pooled_data), "\n")

rm(predictor_data, pixel_ids, stopover_vector, year_vector)
gc()

data_creation_time <- difftime(Sys.time(), start_time, units = "mins")
cat("\nDataset creation completed in", round(data_creation_time, 2), "minutes\n")
cat("Total observations:", nrow(pooled_data), "\n")
cat("Total variables:", ncol(pooled_data), "\n")

#################################
# DATA CLEANING FOR GAM 
#################################

cat("\nCleaning dataset for GAM analysis...\n")
initial_rows <- nrow(pooled_data)

# Remove rows with NA response variable
pooled_data <- pooled_data[!is.na(pooled_data$stopover), ]
cat("Removed", initial_rows - nrow(pooled_data), "rows with NA stopover values\n")

# AUDIT 5: After removing NA stopover
cat("\n=== AUDIT 5: After removing NA stopover ===\n")
audit_log$after_na_stopover <- table(pooled_data$ecoregion)
print_eco_table(pooled_data$ecoregion, "BCR distribution:")
cat("Total observations:", nrow(pooled_data), "\n")

# Get predictor columns
predictor_cols <- names(pooled_data)[!names(pooled_data) %in% c("pixel_id", "stopover", "year")]
cat("Predictor variables:", length(predictor_cols), "\n")

# Remove predictors with all NA values
all_na_predictors <- sapply(pooled_data[predictor_cols], function(x) all(is.na(x)))
if(any(all_na_predictors)) {
  cat("Removing", sum(all_na_predictors), "predictors with all NA values\n")
  pooled_data <- pooled_data[, !names(pooled_data) %in% predictor_cols[all_na_predictors]]
  predictor_cols <- predictor_cols[!all_na_predictors]
}

# Remove predictors with no variation
no_variation <- sapply(pooled_data[predictor_cols], function(x) {
  unique_vals <- unique(x[!is.na(x)])
  length(unique_vals) <= 1
})

if(any(no_variation)) {
  cat("Removing", sum(no_variation), "predictors with no variation\n")
  pooled_data <- pooled_data[, !names(pooled_data) %in% predictor_cols[no_variation]]
  predictor_cols <- predictor_cols[!no_variation]
}

# Handle remaining NA values in predictors
complete_cases <- complete.cases(pooled_data[predictor_cols])
cat("Complete cases:", sum(complete_cases), "of", nrow(pooled_data), "\n")

# AUDIT 6: Before complete.cases filtering
cat("\n=== AUDIT 6: Before complete.cases filtering ===\n")
cat("Will lose:", nrow(pooled_data) - sum(complete_cases), "observations\n")
eco_before <- sort(unique(pooled_data$ecoregion))
eco_after <- sort(unique(pooled_data$ecoregion[complete_cases]))
lost_ecos <- setdiff(eco_before, eco_after)
if(length(lost_ecos) > 0) {
  cat("WARNING: These BCR codes will be completely removed:\n")
  print(bcr_lookup[bcr_lookup$code %in% lost_ecos, ])
}

if(sum(complete_cases) < 1000) {
  cat("WARNING: Very few complete cases. Consider imputation or different approach.\n")
} else {
  pooled_data <- pooled_data[complete_cases, ]
}

# AUDIT 7: After complete.cases
cat("\n=== AUDIT 7: After complete.cases filtering ===\n")
audit_log$after_complete_cases <- table(pooled_data$ecoregion)
print_eco_table(pooled_data$ecoregion, "BCR distribution:")
cat("Total observations:", nrow(pooled_data), "\n")

cat("Final dataset:", nrow(pooled_data), "observations,", length(predictor_cols), "predictors\n")

#################################
# Format data prior to model fitting
#################################

cat("\nStandardizing variables for GAM analysis...\n")

categorical_vars <- c("ecoregion", "NLCD_1995_2020")
spatial_vars <- c("alan_dist_stack", "dcoast")
numeric_vars <- setdiff(predictor_cols, c(categorical_vars, spatial_vars))

cat("Categorical variables:", length(categorical_vars), "\n")
cat("Spatial variables (will center only):", length(spatial_vars), "\n") 
cat("Numeric variables (will standardize):", length(numeric_vars), "\n")

original_stats <- list()

# Standardize numeric variables
for(var in numeric_vars) {
  if(var %in% names(pooled_data)) {
    original_stats[[var]] <- list(
      mean = mean(pooled_data[[var]], na.rm = TRUE),
      sd = sd(pooled_data[[var]], na.rm = TRUE)
    )
    pooled_data[[var]] <- scale(pooled_data[[var]])[,1]
    cat("  Standardized:", var, "\n")
  }
}

# Center spatial variables
for(var in spatial_vars) {
  if(var %in% names(pooled_data)) {
    original_stats[[var]] <- list(
      mean = mean(pooled_data[[var]], na.rm = TRUE),
      sd = sd(pooled_data[[var]], na.rm = TRUE)
    )
    pooled_data[[var]] <- pooled_data[[var]] - mean(pooled_data[[var]], na.rm = TRUE)
    cat("  Centered:", var, "\n")
  }
}

# Convert categorical variables to factors
for(var in categorical_vars) {
  if(var %in% names(pooled_data)) {
    pooled_data[[var]] <- as.factor(pooled_data[[var]])
    cat("  Converted to factor:", var, "with", nlevels(pooled_data[[var]]), "levels\n")
  }
}

# AUDIT 8: After factor conversion
cat("\n=== AUDIT 8: After factor conversion ===\n")
cat("Factor levels (as character):", levels(pooled_data$ecoregion), "\n")
cat("Factor levels (as numeric):", as.numeric(levels(pooled_data$ecoregion)), "\n")
audit_log$factor_table <- table(pooled_data$ecoregion)
print(audit_log$factor_table)
factor_codes <- as.numeric(levels(pooled_data$ecoregion))
cat("\nFactor level to BCR name mapping:\n")
print(bcr_lookup[bcr_lookup$code %in% factor_codes, ])

# Drop unused levels
pooled_data$ecoregion <- droplevels(pooled_data$ecoregion)
pooled_data$NLCD_1995_2020 <- droplevels(pooled_data$NLCD_1995_2020)

# AUDIT 9: After droplevels
cat("\n=== AUDIT 9: After droplevels ===\n")
cat("Final factor levels:", levels(pooled_data$ecoregion), "\n")
print_eco_table(pooled_data$ecoregion, "Final BCR distribution:")

numeric_vars <- setdiff(predictor_cols, c(categorical_vars, spatial_vars))

preprocessing_info <- list(
  original_stats = original_stats,
  categorical_vars = categorical_vars,
  spatial_vars = spatial_vars,
  numeric_vars = numeric_vars,
  sample_size = nrow(pooled_data),
  n_predictors = length(predictor_cols)
)

saveRDS(preprocessing_info, file.path(output_dir, "preprocessing_info.rds"))

# Save audit log
saveRDS(audit_log, file.path(output_dir, "data_audit_log.rds"))
cat("\n=== AUDIT COMPLETE ===\n")
cat("Audit log saved to:", file.path(output_dir, "data_audit_log.rds"), "\n\n")

#################################
# Define Models
#################################

cat("\n", paste(rep("=", 50), collapse=""), "\n")
cat("DEFINING CUSTOM MODELS FOR", toupper(season), "SEASON\n")
cat(paste(rep("=", 50), collapse=""), "\n")

# Check available variables

available_vars <- names(pooled_data)
cat("Available variables:", paste(head(available_vars, 15), collapse=", "), "...\n\n")

# Defining models

# Defining models

if(season == "fall") {
  
  models <- list(
    
    # Model 1: NULL MODEL
    null = as.formula("stopover ~ 1"),
    
    # Model 2: DROUGHT MODEL (5 terms)
    drought = as.formula("stopover ~ 
                        s(SPEI_1995_2020) + 
                        s(SPEI_prev_fall_yr) + 
                        s(SPEI_prev_spring_szn) + 
                        s(SPEI_prev_summer) +
                        s(year, bs='re')"),
    
    # Model 3: CLIMATE MODEL (14 terms)
    climate = as.formula("stopover ~ 
                        s(SPEI_1995_2020) + 
                        s(SPEI_prev_fall_yr) + 
                        s(SPEI_prev_spring_szn) + 
                        s(SPEI_prev_summer) +
                        s(temp_fall_PRISM_stack) + 
                        s(ppt_fall_PRISM_stack) + 
                        s(wetness_fall) +
                        s(fall_uwind) +
                        s(fall_vwind) +
                        ti(SPEI_prev_summer, ppt_fall_PRISM_stack, k=c(4,4)) +
                        ti(SPEI_prev_summer, temp_fall_PRISM_stack, k=c(4,4)) +
                        ti(SPEI_1995_2020, wetness_fall, k=c(4,4)) +
                        ti(ppt_fall_PRISM_stack, temp_fall_PRISM_stack, k=c(4,4)) +
                        s(year, bs='re')"),
    
    # Model 4: HABITAT MODEL (18 terms)
    habitat = as.formula("stopover ~ 
                        s(SPEI_1995_2020) + 
                        s(SPEI_prev_fall_yr) + 
                        s(SPEI_prev_spring_szn) + 
                        s(SPEI_prev_summer) +
                        s(alan_dist_stack) + 
                        s(dcoast) + 
                        s(dwater_perm_stack) + 
                        s(forest_age_stack) +
                        s(msavi_fall) +
                        NLCD_1995_2020 +  
                        s(SPEI_1995_2020, NLCD_1995_2020, bs = 'fs', k=4) +
                        s(SPEI_1995_2020, ecoregion, bs = 'fs', k=4) +
                        ti(SPEI_prev_summer, alan_dist_stack, k=c(4,4)) +
                        ti(SPEI_1995_2020, alan_dist_stack, k=c(4,4)) +
                        ti(SPEI_prev_summer, msavi_fall, k=c(4,4)) +
                        ti(msavi_fall, alan_dist_stack, k=c(4,4)) +
                        s(year, bs='re')"),
    
    # Model 5: FULL MODEL (27 terms)
    full = as.formula("stopover ~ 
                     s(SPEI_1995_2020) +
                     s(SPEI_prev_fall_yr) + 
                     s(SPEI_prev_spring_szn) + 
                     s(SPEI_prev_summer) + 
                     s(temp_fall_PRISM_stack) + 
                     s(ppt_fall_PRISM_stack) + 
                     s(wetness_fall) +
                     s(fall_uwind) +
                     s(fall_vwind) +
                     s(alan_dist_stack) + 
                     s(dcoast) + 
                     s(dwater_perm_stack) + 
                     s(forest_age_stack) +
                     s(msavi_fall) +
                     NLCD_1995_2020 +
                     s(SPEI_1995_2020, NLCD_1995_2020, bs = 'fs', k=4) +
                     s(SPEI_1995_2020, ecoregion, bs = 'fs', k=4) +
                     s(SPEI_prev_summer, ecoregion, bs = 'fs', k=4) +
                     ti(SPEI_prev_summer, ppt_fall_PRISM_stack, k=c(4,4)) +
                     ti(SPEI_prev_summer, temp_fall_PRISM_stack, k=c(4,4)) + 
                     ti(SPEI_1995_2020, wetness_fall, k=c(4,4)) +
                     ti(ppt_fall_PRISM_stack, temp_fall_PRISM_stack, k=c(4,4)) +
                     ti(SPEI_prev_summer, alan_dist_stack, k=c(4,4)) +
                     ti(SPEI_1995_2020, alan_dist_stack, k=c(4,4)) + 
                     ti(SPEI_prev_summer, msavi_fall, k=c(4,4)) +
                     s(year, bs='re')")
  )
  
} else {
  # SPRING SEASON MODELS
  
  models <- list(
    
    # Model 1: NULL MODEL
    null = as.formula("stopover ~ 1"),
    
    # Model 2: DROUGHT MODEL (5 terms)
    drought = as.formula("stopover ~ 
                        s(SPEI_spring_1995_2020) + 
                        s(SPEI_prev_spring_yr) + 
                        s(SPEI_prev_fall_szn) + 
                        s(SPEI_prev_winter) +
                        s(year, bs='re')"),
    
    # Model 3: CLIMATE MODEL (14 terms)
    climate = as.formula("stopover ~ 
                        s(SPEI_spring_1995_2020) + 
                        s(SPEI_prev_spring_yr) + 
                        s(SPEI_prev_fall_szn) + 
                        s(SPEI_prev_winter) + 
                        s(temp_spring_PRISM_stack) + 
                        s(ppt_spring_PRISM_stack) + 
                        s(wetness_spring) +
                        s(uwind_spring) +
                        s(vwind_spring) +
                        ti(SPEI_prev_winter, ppt_spring_PRISM_stack, k=c(4,4)) +
                        ti(SPEI_prev_winter, temp_spring_PRISM_stack, k=c(4,4)) +
                        ti(SPEI_spring_1995_2020, wetness_spring, k=c(4,4)) +
                        ti(ppt_spring_PRISM_stack, temp_spring_PRISM_stack, k=c(4,4)) +
                        s(year, bs='re')"),
    
    # Model 4: HABITAT MODEL (19 terms)
    habitat = as.formula("stopover ~ 
                        s(SPEI_spring_1995_2020) + 
                        s(SPEI_prev_spring_yr) + 
                        s(SPEI_prev_fall_szn) + 
                        s(SPEI_prev_winter) +
                        s(alan_dist_stack) + 
                        s(dcoast) + 
                        s(dwater_perm_stack) + 
                        s(dwater_temp_stack) +
                        s(forest_age_stack) +
                        s(msavi_spring) +
                        NLCD_1995_2020 +
                        s(SPEI_spring_1995_2020, NLCD_1995_2020, bs = 'fs', k=4) +
                        s(SPEI_spring_1995_2020, ecoregion, bs = 'fs', k=4) +
                        ti(SPEI_prev_winter, dwater_temp_stack, k=c(4,4)) +
                        ti(SPEI_spring_1995_2020, dwater_temp_stack, k=c(4,4)) + 
                        ti(SPEI_prev_winter, msavi_spring, k=c(4,4)) +
                        ti(dwater_temp_stack, msavi_spring, k=c(4,4)) +
                        s(year, bs='re')"),
    
    # Model 5: FULL MODEL (28 terms)
    full = as.formula("stopover ~ 
                     s(SPEI_spring_1995_2020) + 
                     s(SPEI_prev_spring_yr) + 
                     s(SPEI_prev_fall_szn) + 
                     s(SPEI_prev_winter) +
                     s(temp_spring_PRISM_stack) + 
                     s(ppt_spring_PRISM_stack) + 
                     s(wetness_spring) +
                     s(uwind_spring) +
                     s(vwind_spring) +
                     s(alan_dist_stack) + 
                     s(dcoast) + 
                     s(dwater_perm_stack) + 
                     s(dwater_temp_stack) +
                     s(forest_age_stack) +
                     s(msavi_spring) +
                     NLCD_1995_2020 +
                     s(SPEI_spring_1995_2020, NLCD_1995_2020, bs = 'fs', k=4) +
                     s(SPEI_spring_1995_2020, ecoregion, bs = 'fs', k=4) +
                     s(SPEI_prev_winter, ecoregion, bs = 'fs', k=4) +
                     ti(SPEI_prev_winter, ppt_spring_PRISM_stack, k=c(4,4)) +
                     ti(SPEI_prev_winter, temp_spring_PRISM_stack, k=c(4,4)) +
                     ti(SPEI_spring_1995_2020, wetness_spring, k=c(4,4)) +
                     ti(ppt_spring_PRISM_stack, temp_spring_PRISM_stack, k=c(4,4)) +
                     ti(SPEI_prev_winter, dwater_temp_stack, k=c(4,4)) +
                     ti(SPEI_spring_1995_2020, dwater_temp_stack, k=c(4,4)) +
                     ti(SPEI_prev_winter, msavi_spring, k=c(4,4)) +
                     s(year, bs='re')")
  )
}
#################################
# SETUP PARALLEL BACKEND (MOVED HERE - AFTER MODELS DEFINED)
#################################

cat("\n", paste(rep("=", 50), collapse=""), "\n")
cat("SETTING UP PARALLEL PROCESSING\n")
cat(paste(rep("=", 50), collapse=""), "\n")

n_cores <- 24  # CHANGE THIS to match your system

# Create cluster - use one worker per model (up to n_cores)
cl <- makeCluster(min(n_cores, length(models)))

# Calculate threads per model
threads_per_model <- max(1, floor(n_cores / length(models)))
cat("Using", length(cl), "workers with", threads_per_model, "threads each\n")

# CRITICAL: Register the cluster as the parallel backend for foreach
registerDoParallel(cl)

# Verify registration
cat("Parallel backend registered:", getDoParName(), "\n")
cat("Number of workers:", getDoParWorkers(), "\n")

# Export necessary objects to ALL workers in the cluster
clusterExport(cl, c("pooled_data", "models", "threads_per_model"), envir = environment())

# Load required packages on ALL workers
clusterEvalQ(cl, {
  library(mgcv)
  library(dplyr)
})

cat("All workers initialized successfully\n")

#################################
# FIT MODELS IN PARALLEL
#################################

cat("\n", paste(rep("=", 50), collapse=""), "\n")
cat("FITTING", length(models), "MODELS IN PARALLEL\n")
cat(paste(rep("=", 50), collapse=""), "\n")

model_fit_start <- Sys.time()

# This loop runs in PARALLEL because we registered the cluster
results <- foreach(
  model_name = names(models),
  .packages = c('mgcv', 'dplyr'),
  .errorhandling = 'pass',  # Continue even if one model fails
  .verbose = FALSE
) %dopar% {  # %dopar% uses the registered parallel backend (cl)
  
  model_start <- Sys.time()
  
  tryCatch({
    
    # Fit model based on complexity
    if(model_name == "null") {
      
      # Null model - use gam()
      fitted_model <- gam(
        formula = models[[model_name]],
        data = pooled_data,
        family = nb(),
        method = "REML"
      )
      
    } else {
      
      fitted_model <- bam(
        formula = models[[model_name]],
        data = pooled_data,
        family = nb(),
        method = "fREML",        
        discrete = TRUE,
	select = TRUE,       
        nthreads = threads_per_model,  # Auto-calculated based on n_cores
        gc.level = 1,           
        control = gam.control(
          trace = FALSE,
          epsilon = 1e-6,
          maxit = 200
        ),
        samfrac = 0.1,        
        chunk.size = 10000      
      )
    }
    
    # Calculate fitting time
    fitting_time <- as.numeric(difftime(Sys.time(), model_start, units = "mins"))
    
    # Extract model statistics
    model_summary <- summary(fitted_model)
    
    # Return results as a list
    list(
      model = fitted_model,
      stats = data.frame(
        model = model_name,
        aic = AIC(fitted_model),
        r_squared = model_summary$r.sq,
        dev_explained = model_summary$dev.expl,
        gcv_score = fitted_model$gcv.ubre,
        edf = sum(fitted_model$edf),
        n_terms = length(fitted_model$coefficients),
        fitting_time = fitting_time,
        stringsAsFactors = FALSE
      ),
      success = TRUE
    )
    
  }, error = function(e) {
    
    # Return error information
    list(
      model = NULL,
      stats = data.frame(
        model = model_name,
        aic = NA,
        r_squared = NA,
        dev_explained = NA,
        gcv_score = NA,
        edf = NA,
        n_terms = NA,
        fitting_time = as.numeric(difftime(Sys.time(), model_start, units = "mins")),
        stringsAsFactors = FALSE
      ),
      success = FALSE,
      error = e$message
    )
  })
}

total_fit_time <- difftime(Sys.time(), model_fit_start, units = "mins")
cat("\nAll models fitted in", round(total_fit_time, 2), "minutes\n")

#################################
# PROCESS RESULTS
#################################

cat("\n", paste(rep("=", 50), collapse=""), "\n")
cat("PROCESSING RESULTS\n")
cat(paste(rep("=", 50), collapse=""), "\n")

fitted_models <- list()
model_results <- data.frame()

for(i in seq_along(results)) {
  model_name <- names(models)[i]
  
  if(results[[i]]$success) {
    fitted_models[[model_name]] <- results[[i]]$model
    model_results <- rbind(model_results, results[[i]]$stats)
    
    cat("✓", toupper(model_name), "- AIC:", round(results[[i]]$stats$aic, 2),
        "R²:", round(results[[i]]$stats$r_squared, 3),
        "Time:", round(results[[i]]$stats$fitting_time, 2), "min\n")
  } else {
    cat("✗", toupper(model_name), "- FAILED:", results[[i]]$error, "\n")
    model_results <- rbind(model_results, results[[i]]$stats)
  }
}

# CRITICAL: Stop the cluster when done
stopCluster(cl)
cat("\nParallel cluster stopped\n")

#################################
# MODEL COMPARISON AND RESULTS
#################################

cat("\n", paste(rep("=", 50), collapse=""), "\n")
cat("MODEL COMPARISON RESULTS\n")
cat(paste(rep("=", 50), collapse=""), "\n")

# Sort by AIC (best model first)

model_results <- model_results[order(model_results$aic, na.last = TRUE), ]

# Calculate AIC differences

model_results$delta_aic <- model_results$aic - min(model_results$aic, na.rm = TRUE)

# Calculate AIC weights

model_results$aic_weight <- exp(-0.5 * model_results$delta_aic) / sum(exp(-0.5 * model_results$delta_aic), na.rm = TRUE)

# Print results table

cat("\nMODEL COMPARISON TABLE (sorted by AIC):\n")
print(model_results[, c("model", "aic", "delta_aic", "aic_weight", "r_squared", "dev_explained", "edf", "fitting_time")])

# Identify best model

best_model <- model_results$model[1]
cat("\nBEST MODEL:", toupper(best_model), "\n")
cat("AIC:", round(model_results$aic[1], 2), "\n")
cat("R²:", round(model_results$r_squared[1], 3), "\n")
cat("Deviance explained:", round(model_results$dev_explained[1], 3), "\n")

#################################
# DETAILED ANALYSIS OF BEST MODEL
#################################

#if(best_model %in% names(fitted_models)) {
 # cat("\n", paste(rep("=", 40), collapse=""), "\n")
 # cat("DETAILED ANALYSIS OF BEST MODEL:", toupper(best_model), "\n")
 # cat(paste(rep("=", 40), collapse=""), "\n")
  
#  best_fitted <- fitted_models[[best_model]]
  
#  cat("\nModel Summary:\n")
#  print(summary(best_fitted))
  
#  cat("\nModel Diagnostics:\n")
#  gam_check_output <- capture.output(gam.check(best_fitted))
#  cat(paste(gam_check_output, collapse="\n"))
#}

#################################
# Saving results
#################################

cat("\nSaving model comparison results...\n")

# Save all fitted models

models_file <- file.path(output_dir, paste0(season, "_all_models1.rds"))
saveRDS(fitted_models, models_file)
cat("All models saved:", models_file, "\n")

# Save comparison results

comparison_file <- file.path(output_dir, paste0(season, "_model_comparison1.csv"))
write.csv(model_results, comparison_file, row.names = FALSE)
cat("Comparison results saved:", comparison_file, "\n")

# Save best model separately

if(best_model %in% names(fitted_models)) {
  best_model_file <- file.path(output_dir, paste0(season, "_best_model1.rds"))
  saveRDS(fitted_models[[best_model]], best_model_file)
  cat("Best model saved:", best_model_file, "\n")
}

# Save dataset

data_file <- file.path(output_dir, paste0(season, "_gam_data.csv"))
write.csv(pooled_data, data_file, row.names = FALSE)
cat("Dataset saved:", data_file, "\n")

total_time <- difftime(Sys.time(), start_time, units = "mins")
cat("\nTotal analysis time:", round(total_time, 2), "minutes\n")
cat("GAM model comparison complete! Results saved in:", output_dir, "\n")
