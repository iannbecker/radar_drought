#################################
#
#   Random Forest
#   Ian Becker
#   9/15/2025
#
################################

library(terra)
library(randomForest)
library(dplyr)
library(doParallel)
library(foreach)

# Set up parallel backend

n_cores <- 8
registerDoParallel(cores = n_cores)

# Set working directory

setwd("/home/ianbecker01/drought_radar")

# Configuration

season <- "spring"  
include_interactions <- TRUE  # Set to FALSE for main effects only
max_interaction_order <- 2    # 2 = pairwise interactions,

# Input directories  

data_dir <- paste0(season, "_raster")  
output_dir <- paste0("random_forest_results/", season)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat(paste(rep("=", 60), collapse=""), "\n")
cat("RANDOM FOREST ANALYSIS -", toupper(season), "SEASON\n")
cat(paste(rep("=", 60), collapse=""), "\n")

# Load raster stacks

cat("Loading raster stacks...\n")

# Stopover data (response variable)

stopover_file <- file.path(data_dir, paste0(season, "_stack.tif"))
if(!file.exists(stopover_file)) {
  # Try alternative naming
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
  name <- gsub(paste0("_", season, "_cropped"), "", name)  # Clean names
  
  cat("Loading predictor:", name, "\n")
  raster_stack <- rast(file)
  predictor_stacks[[name]] <- raster_stack
  cat("  Loaded:", name, "with", nlyr(raster_stack), "layers\n")
}

cat("\nTotal predictor stacks:", length(predictor_stacks), "\n")

# Sample pixels for analysis

cat("\nSampling pixels for model training...\n")
set.seed(123)

# Get valid pixels (non-NA in stopover data)

valid_pixels <- which(!is.na(values(stopover_stack[[1]])))
sample_size <- min(15000, length(valid_pixels))  # Max 25k pixels for memory
sampled_pixels <- sample(valid_pixels, sample_size)

cat("Total valid pixels:", length(valid_pixels), "\n")
cat("Sample size for training:", length(sampled_pixels), "\n")

# Create long-form dataset - EFFICIENT VERSION
cat("\nCreating long-form dataset (efficient approach)...\n")
start_time <- Sys.time()

# First, extract all predictor values upfront
cat("Extracting all predictor values...\n")
predictor_data <- list()
n_years <- nlyr(stopover_stack)

for(pred_name in names(predictor_stacks)) {
  cat("  Extracting values for:", pred_name, "\n")
  pred_stack <- predictor_stacks[[pred_name]]
  
  if(nlyr(pred_stack) == 1) {
    # Static predictor - extract once, replicate across years
    static_vals <- values(pred_stack)[sampled_pixels]
    predictor_data[[pred_name]] <- rep(static_vals, n_years)
    cat("    Static predictor - replicated across", n_years, "years\n")
    
  } else if(nlyr(pred_stack) == n_years) {
    # Time-varying predictor matching stopover years exactly
    time_vals <- values(pred_stack)[sampled_pixels, ]
    predictor_data[[pred_name]] <- as.vector(t(time_vals))  # transpose to get pixel-year order
    cat("    Time-varying predictor -", nlyr(pred_stack), "layers matched\n")
    
  } else {
    # Multiple layers - handle temporal matching
    layer_names <- names(pred_stack)
    
    # Try to extract years from layer names
    layer_years <- regmatches(layer_names, regexpr("\\d{4}", layer_names))
    
    if(length(layer_years) > 0) {
      # Match by year
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
      # Use temporal indexing or create separate predictors
      if(nlyr(pred_stack) >= n_years) {
        # Use first n_years layers
        time_vals <- values(pred_stack[[1:n_years]])[sampled_pixels, ]
        predictor_data[[pred_name]] <- as.vector(t(time_vals))
        cat("    Multi-layer predictor - used first", n_years, "layers\n")
        
      } else {
        # Create separate predictors for each layer
        for(layer_idx in 1:nlyr(pred_stack)) {
          layer_name <- paste0(pred_name, "_L", layer_idx)
          static_vals <- values(pred_stack[[layer_idx]])[sampled_pixels]
          predictor_data[[layer_name]] <- rep(static_vals, n_years)
        }
        cat("    Created", nlyr(pred_stack), "separate predictors\n")
      }
    }
  }
  
  # Clean up immediately
  rm(pred_stack)
  gc()
}

# Extract stopover values (response variable)
cat("Extracting stopover values...\n")
stopover_matrix <- values(stopover_stack)[sampled_pixels, ]
stopover_vector <- as.vector(t(stopover_matrix))

# Clean up stopover stack
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

# Add all predictors
for(pred_name in names(predictor_data)) {
  pooled_data[[pred_name]] <- predictor_data[[pred_name]]
}

# Clean up intermediate objects
rm(predictor_data, pixel_ids, stopover_vector, year_vector)
gc()

data_creation_time <- difftime(Sys.time(), start_time, units = "mins")
cat("\nDataset creation completed in", round(data_creation_time, 2), "minutes\n")
cat("Total observations:", nrow(pooled_data), "\n")
cat("Total variables:", ncol(pooled_data), "\n")

# Data cleaning

cat("\nCleaning dataset...\n")
initial_rows <- nrow(pooled_data)

# Remove rows with NA response variable

pooled_data <- pooled_data[!is.na(pooled_data$stopover), ]
cat("Removed", initial_rows - nrow(pooled_data), "rows with NA stopover values\n")

# Check predictor variables

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

if(sum(complete_cases) < 1000) {
  cat("WARNING: Very few complete cases. Consider imputation or different approach.\n")
} else {
  # Use only complete cases
  pooled_data <- pooled_data[complete_cases, ]
}

cat("Final dataset:", nrow(pooled_data), "observations,", length(predictor_cols), "predictors\n")

# Create interaction terms

if(include_interactions && length(predictor_cols) > 1) {
  cat("\nCreating interaction terms...\n")
  
  # Focus on key predictors for interactions to avoid explosion
  # Prioritize drought, climate, and habitat predictors
  drought_vars <- predictor_cols[grepl("SPEI|drought", predictor_cols, ignore.case = TRUE)]
  climate_vars <- predictor_cols[grepl("ppt|temp|climate|wetness", predictor_cols, ignore.case = TRUE)]
  habitat_vars <- predictor_cols[grepl("NLCD|forest|water|eco|msavi", predictor_cols, ignore.case = TRUE)]
  spatial_vars <- predictor_cols[grepl("alan_dist|dcoast", predictor_cols, ignore.case = TRUE)]
  
  # Always include year in interactions
  key_predictors <- unique(c("year", drought_vars, climate_vars, habitat_vars, spatial_vars))
  key_predictors <- key_predictors[key_predictors %in% predictor_cols]
  
  cat("Creating interactions for", length(key_predictors), "key predictors\n")
  print(key_predictors)
  
  interaction_data <- pooled_data
  n_interactions_created <- 0
  
  if(max_interaction_order >= 2) {
    # Create pairwise interactions
    for(i in 1:(length(key_predictors)-1)) {
      for(j in (i+1):length(key_predictors)) {
        var1 <- key_predictors[i]
        var2 <- key_predictors[j]
        
        # Skip if either variable is categorical with too many levels
        if(is.factor(interaction_data[[var1]]) && nlevels(interaction_data[[var1]]) > 10) next
        if(is.factor(interaction_data[[var2]]) && nlevels(interaction_data[[var2]]) > 10) next
        
        interaction_name <- paste0(var1, "_X_", var2)
        
        # Create interaction term
        if(is.numeric(interaction_data[[var1]]) && is.numeric(interaction_data[[var2]])) {
          # Numeric * Numeric
          interaction_data[[interaction_name]] <- interaction_data[[var1]] * interaction_data[[var2]]
          n_interactions_created <- n_interactions_created + 1
          
        } else if(is.factor(interaction_data[[var1]]) || is.factor(interaction_data[[var2]])) {
          # Categorical interactions - create dummy interactions
          if(is.factor(interaction_data[[var1]])) {
            # Factor * Numeric
            for(level in levels(interaction_data[[var1]])[1:min(5, nlevels(interaction_data[[var1]]))]) {
              dummy_name <- paste0(interaction_name, "_", level)
              dummy_var <- as.numeric(interaction_data[[var1]] == level)
              interaction_data[[dummy_name]] <- dummy_var * interaction_data[[var2]]
              n_interactions_created <- n_interactions_created + 1
            }
          } else {
            # Numeric * Factor
            for(level in levels(interaction_data[[var2]])[1:min(5, nlevels(interaction_data[[var2]]))]) {
              dummy_name <- paste0(interaction_name, "_", level)
              dummy_var <- as.numeric(interaction_data[[var2]] == level)
              interaction_data[[dummy_name]] <- interaction_data[[var1]] * dummy_var
              n_interactions_created <- n_interactions_created + 1
            }
          }
        }
        
        # Prevent memory explosion
        if(n_interactions_created > 500) {
          cat("Stopping interaction creation at 500 terms to prevent memory issues\n")
          break
        }
      }
      if(n_interactions_created > 500) break
    }
  }
  
  cat("Created", n_interactions_created, "interaction terms\n")
  
  # Update dataset
  pooled_data <- interaction_data
  
  # Update predictor columns list
  new_predictors <- setdiff(names(pooled_data), c("pixel_id", "stopover", "year"))
  cat("Total predictors including interactions:", length(new_predictors), "\n")
  
  # Clean interaction terms - remove those with no variation
  interaction_cols <- names(pooled_data)[grepl("_X_", names(pooled_data))]
  if(length(interaction_cols) > 0) {
    no_var_interactions <- sapply(pooled_data[interaction_cols], function(x) {
      unique_vals <- unique(x[!is.na(x)])
      length(unique_vals) <= 1
    })
    
    if(any(no_var_interactions)) {
      cat("Removing", sum(no_var_interactions), "interaction terms with no variation\n")
      pooled_data <- pooled_data[, !names(pooled_data) %in% interaction_cols[no_var_interactions]]
    }
  }
}

# Data summary

cat("\nResponse variable summary:\n")
cat("Stopover - Mean:", round(mean(pooled_data$stopover), 3), 
    "SD:", round(sd(pooled_data$stopover), 3), "\n")
cat("Year range:", range(pooled_data$year), "\n")

# Train Random Forest model

cat("\n", paste(rep("=", 40), collapse=""), "\n")
cat("TRAINING RANDOM FOREST MODEL\n")
cat(paste(rep("=", 40), collapse=""), "\n")

rf_start_time <- Sys.time()

# Prepare data for RF - get final predictor list

final_predictors <- setdiff(names(pooled_data), c("pixel_id", "stopover"))
rf_data <- pooled_data[, c("stopover", final_predictors)]

# Set RF parameters

ntree <- 1000

# Increase mtry for interactions to sample more variables per split

mtry <- max(1, floor(sqrt(ncol(rf_data) - 1) * 1.5))  # -1 to exclude response
mtry <- min(mtry, ncol(rf_data) - 1)  # Cap at max possible

cat("Training Random Forest with:\n")
cat("- Observations:", nrow(rf_data), "\n")
cat("- Predictors:", ncol(rf_data) - 1, "\n") # -1 for response
cat("- Trees:", ntree, "\n")
cat("- Variables per split:", mtry, "\n")

# Train model

rf_model <- foreach(ntree = rep(ceiling(1000/n_cores), n_cores), 
                    .combine = randomForest::combine,
                    .multicombine = TRUE,
                    .packages = "randomForest") %dopar% {
                      randomForest(stopover ~ ., 
                                   data = rf_data,
                                   ntree = ntree,
                                   mtry = mtry,
                                   importance = TRUE,
                                   na.action = na.omit)
                    }

rf_training_time <- difftime(Sys.time(), rf_start_time, units = "mins")

# Model results

# printing results

cat("\n", paste(rep("=", 40), collapse=""), "\n")
cat("RANDOM FOREST RESULTS\n")
cat(paste(rep("=", 40), collapse=""), "\n")

cat("Training time:", round(rf_training_time, 2), "minutes\n")

# Variable importance with interaction analysis

cat("\nVariable Importance Analysis:\n")
importance_df <- data.frame(
  variable = rownames(importance(rf_model)),
  inc_mse = importance(rf_model)[, "%IncMSE"],
  inc_node_purity = importance(rf_model)[, "IncNodePurity"]
)

importance_df <- importance_df[order(importance_df$inc_mse, decreasing = TRUE), ]

# Separate main effects from interactions

main_effects <- importance_df[!grepl("_X_", importance_df$variable), ]
interactions <- importance_df[grepl("_X_", importance_df$variable), ]

cat("\nTop 15 Main Effects:\n")
print(head(main_effects, 15))

if(nrow(interactions) > 0) {
  cat("\nTop 10 Interaction Effects:\n")
  print(head(interactions, 10))
  
  cat("\nInteraction Summary:\n")
  cat("Total interactions tested:", nrow(interactions), "\n")
  cat("Interactions in top 50 variables:", sum(head(importance_df, 50)$variable %in% interactions$variable), "\n")
}

# Variable selection recommendations

cat("\nVariable Selection Recommendations:\n")

# Top variables by importance

top_20_vars <- head(importance_df, 20)$variable
drought_in_top20 <- sum(grepl("SPEI|drought", top_20_vars, ignore.case = TRUE))
climate_in_top20 <- sum(grepl("ppt|temp|wind", top_20_vars, ignore.case = TRUE))
habitat_in_top20 <- sum(grepl("NLCD|forest|water|eco", top_20_vars, ignore.case = TRUE))
interactions_in_top20 <- sum(grepl("_X_", top_20_vars))

cat("In top 20 variables:\n")
cat("- Drought variables:", drought_in_top20, "\n")
cat("- Climate variables:", climate_in_top20, "\n") 
cat("- Habitat variables:", habitat_in_top20, "\n")
cat("- Interaction terms:", interactions_in_top20, "\n")

# Year effect analysis

if("year" %in% importance_df$variable) {
  year_rank <- which(importance_df$variable == "year")
  year_importance <- importance_df$inc_mse[year_rank]
  cat("\nTemporal trend analysis:\n")
  cat("Year importance ranking:", year_rank, "of", nrow(importance_df), "\n")
  cat("Year importance score:", round(year_importance, 3), "\n")
  
  if(year_rank <= 10) {
    cat("-> Strong temporal trend detected - consider time-varying coefficients in GAMs\n")
  } else if(year_rank <= nrow(importance_df) * 0.3) {
    cat("-> Moderate temporal trend - include year as covariate\n")
  } else {
    cat("-> Weak temporal trend - may not need explicit year effects\n")
  }
}

# Save results

cat("\nSaving results...\n")

# Save model

model_file <- file.path(output_dir, "random_forest_model.rds")
saveRDS(rf_model, model_file)
cat("Model saved:", model_file, "\n")

# Save dataset

data_file <- file.path(output_dir, "analysis_dataset.csv")
write.csv(pooled_data, data_file, row.names = FALSE)
cat("Dataset saved:", data_file, "\n")

# Save importance results

importance_file <- file.path(output_dir, "variable_importance.csv")
write.csv(importance_df, importance_file, row.names = FALSE)
cat("Variable importance saved:", importance_file, "\n")

# Save summary

summary_data <- data.frame(
  metric = c("n_observations", "n_predictors", "ntree", "mtry", "oob_mse", "r_squared", "training_time_mins"),
  value = c(nrow(rf_data), ncol(rf_data)-1, ntree, mtry, oob_mse, r_squared, as.numeric(rf_training_time))
)

summary_file <- file.path(output_dir, "model_summary.csv")
write.csv(summary_data, summary_file, row.names = FALSE)
cat("Summary saved:", summary_file, "\n")

total_time <- difftime(Sys.time(), start_time, units = "mins")
cat("\nTotal analysis time:", round(total_time, 2), "minutes\n")
cat("Analysis complete! Results saved in:", output_dir, "\n")

stopImplicitCluster()
