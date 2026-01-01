##############################
#
#   Random Forest - Fall & Spring 
#   Ian Becker
#   December 2025
#
##############################

# This script is used to run random forest for our predictors
# in order to select interaction terms for our GAMs

library(randomForest)
library(doParallel)
library(foreach)

setwd("PATH HERE")

#################################
# PREP
#################################

cat("\nPREPARING DATA FOR RANDOM FOREST\n")

# Select season (change between fall and spring as needed)

season <- "fall"  

# Input/output paths

input_data <- paste0(season, "_rf_data.csv")  # Pre-processed dataset
output_dir <- paste0("random_forest_results/", season)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Random Forest parameters

ntree <- 1000
n_cores <- 8

# Setup parallel processing for later

registerDoParallel(cores = n_cores)

# Loading data

cat("LOADING DATA FOR", toupper(season), "\n")
full_data <- read.csv(input_data, stringsAsFactors = FALSE)

# Inspect data

cat("DATASET LOADED:\n")
cat("TOTAL OBSERVATIONS:", nrow(full_data), "\n")
cat("VARIABLES:", ncol(full_data), "\n")

#################################
# SAMPLE DATA (OPTIONAL)
#################################

sample_size <- NULL  # Set to NULL to use all data

if(!is.null(sample_size) && sample_size < nrow(full_data)) {
  cat("\nSampling data for Random Forest training...\n")
  set.seed(123)
  sample_indices <- sample(1:nrow(full_data), sample_size)
  pooled_data <- full_data[sample_indices, ]
  cat("  Sampled", nrow(pooled_data), "observations from", nrow(full_data), "total\n")
} else {
  pooled_data <- full_data
  cat("  Using all", nrow(pooled_data), "observations\n")
}

# Get predictor columns

predictor_cols <- names(pooled_data)[!names(pooled_data) %in% c("pixel_id", "stopover", "year")]
cat("  Predictor variables:", length(predictor_cols), "\n")

# Convert categorical variables to factors

categorical_vars <- c("ecoregion", "NLCD_1995_2020")
for(var in categorical_vars) {
  if(var %in% names(pooled_data)) {
    pooled_data[[var]] <- as.factor(pooled_data[[var]])
    cat("Converted", var, "to factor with", nlevels(pooled_data[[var]]), "levels\n")
  }
}

#################################
# CREATE INTERACTION TERMS
#################################

cat("\nCreating interaction terms...\n")

# Identify key predictors for interactions

drought_vars <- grep("SPEI|drought", predictor_cols, value = TRUE, ignore.case = TRUE)
climate_vars <- grep("temp|ppt|wind|wetness", predictor_cols, value = TRUE, ignore.case = TRUE)
habitat_vars <- grep("alan|forest|water|msavi|dcoast", predictor_cols, value = TRUE, ignore.case = TRUE)

key_predictors <- unique(c(drought_vars, climate_vars, habitat_vars))
cat("Creating interactions among", length(key_predictors), "key predictors\n")

interaction_data <- pooled_data
n_interactions_created <- 0

# Pairwise interactions

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
    
    # Prevent too many interactions 
    if(n_interactions_created > 500) {
      cat("Stopping interaction creation at 500 terms to prevent memory issues\n")
      break
    }
  }
  if(n_interactions_created > 500) break
}

cat("Created", n_interactions_created, "interaction terms\n")

# Update dataset

pooled_data <- interaction_data

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

cat("Total predictors including interactions:", 
    length(setdiff(names(pooled_data), c("pixel_id", "stopover", "year"))), "\n")

#################################
# TRAIN RANDOM FOREST MODEL
#################################

cat("TRAINING RANDOM FOREST MODEL\n")

# track time 

rf_start_time <- Sys.time()

# Prepare data for RF

final_predictors <- setdiff(names(pooled_data), c("pixel_id", "stopover"))
rf_data <- pooled_data[, c("stopover", final_predictors)]

# Calculate mtry (variables per split)

mtry <- max(1, floor(sqrt(ncol(rf_data) - 1) * 1.5))
mtry <- min(mtry, ncol(rf_data) - 1)

# Data check before running 

cat("Training Random Forest with:\n")
cat("  Observations:", nrow(rf_data), "\n")
cat("  Predictors:", ncol(rf_data) - 1, "\n")
cat("  Trees:", ntree, "\n")
cat("  Variables per split (mtry):", mtry, "\n")

# Train model in parallel (for cluster work)

rf_model <- foreach(ntree = rep(ceiling(ntree/n_cores), n_cores), 
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

#################################
# MODEL RESULTS
#################################

cat("RANDOM FOREST RESULTS\n")

# Check time to train

cat("Training time:", round(rf_training_time, 2), "minutes\n")

# Extract variable importance

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
  cat("  Total interactions tested:", nrow(interactions), "\n")
  cat("  Interactions in top 50 variables:", 
      sum(head(importance_df, 50)$variable %in% interactions$variable), "\n")
}

#################################
# SAVE RESULTS
#################################

cat("\nSaving results...\n")

# Save model

model_file <- file.path(output_dir, paste0(season, "_rf_model.rds"))
saveRDS(rf_model, model_file)
cat("  Model saved:", model_file, "\n")

# Save variable importance

importance_file <- file.path(output_dir, paste0(season, "_variable_importance.csv"))
write.csv(importance_df, importance_file, row.names = FALSE)
cat("  Variable importance saved:", importance_file, "\n")

# Save model summary

summary_data <- data.frame(
  metric = c("n_observations", "n_predictors", "ntree", "mtry", "training_time_mins"),
  value = c(nrow(rf_data), ncol(rf_data)-1, ntree, mtry, as.numeric(rf_training_time))
)

summary_file <- file.path(output_dir, paste0(season, "_model_summary.csv"))
write.csv(summary_data, summary_file, row.names = FALSE)
cat("  Summary saved:", summary_file, "\n")

stopImplicitCluster()

cat("Random Forest analysis complete! Results saved in:", output_dir, "\n")