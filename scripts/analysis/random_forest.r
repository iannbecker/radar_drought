##############################
#
#   Random Forest 
#   Ian Becker
#   May 2026
#
##############################

# This script is used to run the covariate screening random forest analysis.
# Originally, this script was developed to be run on a cluster

library(randomForest)
library(dplyr)
library(doParallel)
library(foreach)

setwd("PATH HERE")

# SCRIPT PREP ------------------------------------

season   <- "fall"
n_cores  <- 8
ntree    <- 1000
mtry     <- 2        # 2 variables per split 

set.seed(123) # seed for reproducibility

source_dir <- "PATH HERE"
output_dir <- "PATH HERE"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# This is for tracking purposes

total_start <- Sys.time()

cat(paste(rep("=", 60), collapse=""), "\n")
cat("RANDOM FOREST -", toupper(season), "SEASON\n")
cat("Spatial subsampling at 4.7km (matches SPEI native resolution)\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# LOAD IN DATA ------------------------------------

# Look for dataset

source_file <- file.path(source_dir, paste0(season, "_source_data.rds"))
lookup_file <- file.path(source_dir, paste0(season, "_pixel_lookup.rds"))

# Error check

if(!file.exists(source_file)) stop("Source data not found: ", source_file,
                                   "\nRun bootstrap_test_single.R first.")
if(!file.exists(lookup_file)) stop("Pixel lookup not found: ", lookup_file)

# Load in data

pooled_data_full <- readRDS(source_file)
pixel_lookup     <- readRDS(lookup_file)

cat("Source data loaded:", nrow(pooled_data_full), "observations\n")
cat("Pixel lookup:", nrow(pixel_lookup), "pixels\n")
cat("Radar IDs:", paste(sort(unique(pixel_lookup$radar_id)), collapse=", "), "\n\n")

# CLEAN DATA ------------------------------------

cat("Cleaning dataset...\n")

# Clean NAs

pooled_data_full <- pooled_data_full[!is.na(pooled_data_full$stopover), ]

# Filter for predictors only 

predictor_cols <- names(pooled_data_full)[!names(pooled_data_full) %in%
                                            c("pixel_id", "stopover", "year")]

# Remove all-NA predictors

all_na <- sapply(pooled_data_full[predictor_cols], function(x) all(is.na(x)))
if(any(all_na)) {
  pooled_data_full <- pooled_data_full[, !names(pooled_data_full) %in%
                                         predictor_cols[all_na]]
  predictor_cols   <- predictor_cols[!all_na]
  cat("Removed", sum(all_na), "all-NA predictors\n")
}

# Remove zero-variance predictors

no_var <- sapply(pooled_data_full[predictor_cols], function(x) {
  length(unique(x[!is.na(x)])) <= 1
})
if(any(no_var)) {
  pooled_data_full <- pooled_data_full[, !names(pooled_data_full) %in%
                                         predictor_cols[no_var]]
  predictor_cols   <- predictor_cols[!no_var]
  cat("Removed", sum(no_var), "zero-variance predictors\n")
}

# Complete cases only

complete_cases   <- complete.cases(pooled_data_full[predictor_cols])
pooled_data_full <- pooled_data_full[complete_cases, ]

cat("Clean dataset:", nrow(pooled_data_full), "observations,",
    length(predictor_cols), "predictors\n")

# Sampling only draws from pixels with complete data

pixel_lookup <- pixel_lookup[
  pixel_lookup$pixel_id %in% unique(pooled_data_full$pixel_id), ]

cat("Updated pixel lookup:", nrow(pixel_lookup), "valid pixels\n\n")

# SPEI RESOLUTION SUBSAMPLING ------------------------------------

# Sample one pixel per 4.7km coarse cell per radar
# Matches native SPEI resolution 

cat("Applying 4.7km spatial subsampling...\n")

sampled_pixel_ids <- unlist(lapply(unique(pixel_lookup$radar_id), function(rid) {
  rid_rows     <- pixel_lookup[pixel_lookup$radar_id == rid, ]
  cell_samples <- tapply(rid_rows$pixel_id, rid_rows$coarse_cell_id,
                         function(px) sample(px, 1))
  unlist(cell_samples)
}))

pooled_data <- pooled_data_full[pooled_data_full$pixel_id %in% sampled_pixel_ids, ]

cat("Sampled pixels:", length(sampled_pixel_ids), "\n")
cat("Observations:", nrow(pooled_data), "\n\n")

rm(pooled_data_full)
gc()

# PREP AND TRAIN RANDOM FOREST MODEL ------------------------------------

# Exclude pixel_id from predictors — keep year as a predictor

final_predictors <- setdiff(names(pooled_data), c("pixel_id", "stopover"))
rf_data          <- pooled_data[, c("stopover", final_predictors)]

cat("RF dataset:\n")
cat("  Observations:", nrow(rf_data), "\n")
cat("  Predictors:", ncol(rf_data) - 1, "\n")
cat("  Trees:", ntree, "\n")
cat("  Variables per split (mtry):", mtry, "\n\n")

# Parallel across cores — combine trees at end

cat(paste(rep("=", 50), collapse=""), "\n")
cat("TRAINING RANDOM FOREST\n")
cat(paste(rep("=", 50), collapse=""), "\n\n")

# Set up parallel backend

registerDoParallel(cores = n_cores)
rf_start <- Sys.time()

trees_per_core <- ceiling(ntree / n_cores)

# Run random forest in parallel

rf_model <- foreach(
  i             = 1:n_cores,
  .combine      = randomForest::combine,
  .multicombine = TRUE,
  .packages     = "randomForest"
) %dopar% {
  randomForest(
    stopover ~ .,
    data       = rf_data,
    ntree      = trees_per_core,
    mtry       = mtry,
    importance = TRUE,
    na.action  = na.omit
  )
}

stopImplicitCluster()

rf_time <- difftime(Sys.time(), rf_start, units = "mins")
cat("Training complete in", round(rf_time, 2), "minutes\n\n")

# VARIABLE IMPORTANCE ------------------------------------

cat(paste(rep("=", 50), collapse=""), "\n")
cat("VARIABLE IMPORTANCE\n")
cat(paste(rep("=", 50), collapse=""), "\n\n")

importance_df <- data.frame(
  variable        = rownames(importance(rf_model)),
  inc_mse         = importance(rf_model)[, "%IncMSE"],
  inc_node_purity = importance(rf_model)[, "IncNodePurity"],
  stringsAsFactors = FALSE
)

importance_df <- importance_df[order(importance_df$inc_mse, decreasing = TRUE), ]

cat("Top 20 variables by % increase in MSE:\n")
print(head(importance_df, 20))

# SAVE RESULTS ------------------------------------

cat("\nSaving results...\n")

saveRDS(rf_model,
        file.path(output_dir, paste0(season, "_rf_model.rds")))

write.csv(importance_df,
          file.path(output_dir, paste0(season, "_variable_importance.csv")),
          row.names = FALSE)

summary_df <- data.frame(
  metric = c("n_observations", "n_predictors", "ntree", "mtry",
             "oob_mse", "oob_r_squared", "training_time_mins"),
  value  = c(nrow(rf_data), ncol(rf_data) - 1, ntree, mtry,
             oob_mse, r_squared, as.numeric(rf_time))
)

write.csv(summary_df,
          file.path(output_dir, paste0(season, "_rf_summary.csv")),
          row.names = FALSE)

total_time <- difftime(Sys.time(), total_start, units = "mins")

cat(paste(rep("=", 60), collapse=""), "\n")
cat("COMPLETE\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("Total time:", round(total_time, 2), "minutes\n")
cat("Results saved to:", output_dir, "\n")