##############################
#
#   GAM Model Comparison
#   Ian Becker
#
##############################

# This script fits 9 competing GAMs for fall or spring migration.
# Originally this script was built to run on a cluster

library(mgcv)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)

setwd("PATH HERE")

# SCRIPT PREP ------------------------------------

season      <- "fall"   # Change to "spring" for spring job
sample_frac <- 0.30     # Proportion of pixels sampled per radar footprint
n_cores     <- 45       # Matches --cpus-per-task in batch script

# Setup paths 

source_dir <- paste0("PATH", season)
output_dir <- paste0("PATH", season)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Tracking 

total_start <- Sys.time()

cat(paste(rep("=", 60), collapse=""), "\n")
cat("GAM MODEL COMPARISON -", toupper(season), "SEASON\n")
cat("Sample fraction:", sample_frac * 100, "%\n")
cat(paste(rep("=", 60), collapse=""), "\n\n")

# LOAD IN DATA ------------------------------------

# Look up data

source_file <- file.path(source_dir, paste0(season, "_source_data.rds"))
lookup_file <- file.path(source_dir, paste0(season, "_pixel_lookup.rds"))

# Error check

if(!file.exists(source_file)) stop("Source data not found: ", source_file,
                                   "\nRun bootstrap_test_single.R first.")
if(!file.exists(lookup_file)) stop("Pixel lookup not found: ", lookup_file)

# Read in data

cat("Loading source data...\n")
pooled_data_full <- readRDS(source_file)
pixel_lookup     <- readRDS(lookup_file)

cat("Source data loaded:", nrow(pooled_data_full), "observations\n")
cat("Pixel lookup:", nrow(pixel_lookup), "pixels\n")
cat("Radar IDs:", paste(sort(unique(pixel_lookup$radar_id)), collapse=", "), "\n\n")

# STRATIFIED SAMPLING BY RADAR ------------------------------------

# Samples a fixed proportion of pixels from each radar footprint
# ensuring equal representation across all available radars.
# Year-level availability is handled downstream via NA stopover removal.

cat("Stratified sampling", sample_frac * 100, "% of pixels per radar...\n")

# Set for reproducibility

set.seed(123)

# Get unique radar IDs, excluding any NA values

radar_ids <- sort(unique(pixel_lookup$radar_id))
radar_ids <- radar_ids[!is.na(radar_ids)]

# Sample pixels within each radar footprint

sampled_pixel_ids <- unlist(lapply(radar_ids, function(rid) {
  px       <- pixel_lookup$pixel_id[pixel_lookup$radar_id == rid]
  n_sample <- max(1, round(length(px) * sample_frac))
  sample(px, min(n_sample, length(px)))
}))

# Subset the full dataset to only include sampled pixels

pooled_data <- pooled_data_full[pooled_data_full$pixel_id %in% sampled_pixel_ids, ]

cat("Sampled pixels:", length(sampled_pixel_ids), "\n")
cat("Observations:", nrow(pooled_data), "\n")
cat("Using", round(length(sampled_pixel_ids) /
                     nrow(pixel_lookup) * 100, 1), "% of valid pixels\n\n")

rm(pooled_data_full)
gc()

# CORRELATION CHECK ------------------------------------

cat("Running predictor correlation check...\n")

# Exclude non-numeric predictor columns

corr_vars <- predictor_cols[!predictor_cols %in%
                              c("ecoregion", "NLCD_1995_2020", "radar_id", "year")]

# Spearman correlation

corr_matrix <- cor(pooled_data[, corr_vars],
                   use = "pairwise.complete.obs", method = "spearman")

# Identify highly correlated pairs

corr_upper <- corr_matrix
corr_upper[lower.tri(corr_upper, diag = TRUE)] <- NA

high_corr <- which(abs(corr_upper) >= 0.7, arr.ind = TRUE)

if(nrow(high_corr) > 0) {
  cat("Highly correlated pairs (|r| >= 0.7):\n")
  high_corr_df <- data.frame(
    var1       = rownames(corr_upper)[high_corr[, 1]],
    var2       = colnames(corr_upper)[high_corr[, 2]],
    spearman_r = corr_upper[high_corr]
  )
  high_corr_df <- high_corr_df[order(abs(high_corr_df$spearman_r),
                                     decreasing = TRUE), ]
  print(high_corr_df)
} else {
  cat("No pairs exceed correlation threshold of 0.7\n")
}

# Save correlation data for reference 

write.csv(as.data.frame(corr_matrix),
          file.path(output_dir, paste0(season, "_predictor_correlations.csv")))

# FORMAT DATA ------------------------------------

cat("\nStandardizing variables...\n")

categorical_vars <- c("ecoregion", "NLCD_1995_2020", "radar_id")
spatial_vars     <- c("alan_dist_stack", "dcoast")
numeric_vars     <- setdiff(predictor_cols, c(categorical_vars, spatial_vars))

original_stats <- list()

# Standardize numeric predictors 

for(var in numeric_vars) {
  if(var %in% names(pooled_data)) {
    original_stats[[var]] <- list(
      mean = mean(pooled_data[[var]], na.rm = TRUE),
      sd   = sd(pooled_data[[var]], na.rm = TRUE)
    )
    pooled_data[[var]] <- scale(pooled_data[[var]])[, 1]
  }
}

# Center spatial predictors 

for(var in spatial_vars) {
  if(var %in% names(pooled_data)) {
    original_stats[[var]] <- list(
      mean = mean(pooled_data[[var]], na.rm = TRUE),
      sd   = sd(pooled_data[[var]], na.rm = TRUE)
    )
    pooled_data[[var]] <- pooled_data[[var]] -
      mean(pooled_data[[var]], na.rm = TRUE)
  }
}

# Convert categorical predictors to factors and drop unused levels

for(var in categorical_vars) {
  if(var %in% names(pooled_data)) {
    pooled_data[[var]] <- as.factor(pooled_data[[var]])
  }
}

# Drop levels from factors

pooled_data$ecoregion      <- droplevels(pooled_data$ecoregion)
pooled_data$NLCD_1995_2020 <- droplevels(pooled_data$NLCD_1995_2020)
pooled_data$radar_id       <- droplevels(pooled_data$radar_id)

cat("Factor levels — ecoregion:", levels(pooled_data$ecoregion), "\n")
cat("Factor levels — radar_id:", levels(pooled_data$radar_id), "\n\n")

# Save original stats for later

saveRDS(list(
  original_stats   = original_stats,
  categorical_vars = categorical_vars,
  spatial_vars     = spatial_vars,
  numeric_vars     = numeric_vars,
  sample_size      = nrow(pooled_data),
  n_predictors     = length(predictor_cols)
), file.path(output_dir, "preprocessing_info.rds"))

# DEFINE MODELS ------------------------------------

cat(paste(rep("=", 50), collapse=""), "\n")
cat("DEFINING MODELS FOR", toupper(season), "SEASON\n")
cat(paste(rep("=", 50), collapse=""), "\n\n")

if(season == "fall") {
  
  models <- list(
    
    # 1: NULL
    
    null = as.formula("stopover ~ 1"),
    
    # 2: DROUGHT — additive only
    
    drought = as.formula("stopover ~
                         s(SPEI_1995_2020) +
                         s(SPEI_prev_fall_yr) +
                         s(SPEI_prev_spring_szn) +
                         s(SPEI_prev_summer) +
                         s(year, bs='re') +
                         s(radar_id, bs='re')"),
    
    # 3: CLIMATE — additive only
    
    climate_additive = as.formula("stopover ~
                         s(SPEI_1995_2020) +
                         s(SPEI_prev_fall_yr) +
                         s(SPEI_prev_spring_szn) +
                         s(SPEI_prev_summer) +
                         s(temp_fall_PRISM_stack) +
                         s(ppt_fall_PRISM_stack) +
                         s(wetness_fall) +
                         s(fall_uwind, k=20) +
                         s(fall_vwind, k=20) +
                         s(year, bs='re') +
                         s(radar_id, bs='re')"),
    
    # 4: CLIMATE — with interactions
    
    climate_interaction = as.formula("stopover ~
                         s(SPEI_1995_2020) +
                         s(SPEI_prev_fall_yr) +
                         s(SPEI_prev_spring_szn) +
                         s(SPEI_prev_summer) +
                         s(temp_fall_PRISM_stack) +
                         s(ppt_fall_PRISM_stack) +
                         s(wetness_fall) +
                         s(fall_uwind, k=20) +
                         s(fall_vwind, k=20) +
                         ti(SPEI_1995_2020, SPEI_prev_summer, k=c(9,9)) +
                         ti(SPEI_1995_2020, temp_fall_PRISM_stack, k=c(9,9)) +
                         ti(fall_uwind, fall_vwind, k=c(9,9)) +
                         s(year, bs='re') +
                         s(radar_id, bs='re')"),
    
    # 5: HABITAT — additive only
    
    habitat_additive = as.formula("stopover ~
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
                         s(year, bs='re') +
                         s(radar_id, bs='re')"),
    
    # 6: HABITAT — with interactions
    
    habitat_interaction = as.formula("stopover ~
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
                         s(SPEI_1995_2020, NLCD_1995_2020, bs='fs', k=4) +
                         s(SPEI_1995_2020, ecoregion, bs='fs', k=4) +
                         s(SPEI_prev_summer, ecoregion, bs='fs', k=4) +
                         ti(SPEI_1995_2020, SPEI_prev_summer, k=c(9,9)) +
                         ti(SPEI_1995_2020, temp_fall_PRISM_stack, k=c(9,9)) +
                         ti(SPEI_prev_summer, alan_dist_stack, k=c(9,9)) +
                         ti(SPEI_1995_2020, alan_dist_stack, k=c(9,9)) +
                         ti(SPEI_prev_summer, msavi_fall, k=c(9,9)) +
                         ti(msavi_fall, alan_dist_stack, k=c(9,9)) +
                         s(year, bs='re') +
                         s(radar_id, bs='re')"),
    
    # 7: DROUGHT INTERACTION — drought terms + concurrent x previous season
    
    drought_interaction = as.formula("stopover ~
                         s(SPEI_1995_2020) +
                         s(SPEI_prev_fall_yr) +
                         s(SPEI_prev_spring_szn) +
                         s(SPEI_prev_summer) +
                         ti(SPEI_1995_2020, SPEI_prev_summer, k=c(9,9)) +
                         s(year, bs='re') +
                         s(radar_id, bs='re')"),
    
    # 8: FULL — additive only
    
    full_additive = as.formula("stopover ~
                         s(SPEI_1995_2020) +
                         s(SPEI_prev_fall_yr) +
                         s(SPEI_prev_spring_szn) +
                         s(SPEI_prev_summer) +
                         s(temp_fall_PRISM_stack) +
                         s(ppt_fall_PRISM_stack) +
                         s(wetness_fall) +
                         s(fall_uwind, k=20) +
                         s(fall_vwind, k=20) +
                         s(alan_dist_stack) +
                         s(dcoast) +
                         s(dwater_perm_stack) +
                         s(forest_age_stack) +
                         s(msavi_fall) +
                         NLCD_1995_2020 +
                         s(year, bs='re') +
                         s(radar_id, bs='re')"),
    
    # 9: FULL — with all interactions
    
    full_interaction = as.formula("stopover ~
                         s(SPEI_1995_2020) +
                         s(SPEI_prev_fall_yr) +
                         s(SPEI_prev_spring_szn) +
                         s(SPEI_prev_summer) +
                         s(temp_fall_PRISM_stack) +
                         s(ppt_fall_PRISM_stack) +
                         s(wetness_fall) +
                         s(fall_uwind, k=20) +
                         s(fall_vwind, k=20) +
                         s(alan_dist_stack) +
                         s(dcoast) +
                         s(dwater_perm_stack) +
                         s(forest_age_stack) +
                         s(msavi_fall) +
                         NLCD_1995_2020 +
                         s(SPEI_1995_2020, NLCD_1995_2020, bs='fs', k=4) +
                         s(SPEI_1995_2020, ecoregion, bs='fs', k=4) +
                         s(SPEI_prev_summer, ecoregion, bs='fs', k=4) +
                         ti(SPEI_1995_2020, SPEI_prev_summer, k=c(9,9)) +
                         ti(SPEI_1995_2020, temp_fall_PRISM_stack, k=c(9,9)) +
                         ti(ppt_fall_PRISM_stack, msavi_fall, k=c(9,9)) +
                         ti(SPEI_prev_summer, alan_dist_stack, k=c(9,9)) +
                         ti(SPEI_1995_2020, alan_dist_stack, k=c(9,9)) +
                         ti(SPEI_prev_summer, msavi_fall, k=c(9,9)) +
                         ti(fall_uwind, fall_vwind, k=c(9,9)) +
                         s(year, bs='re') +
                         s(radar_id, bs='re')")
  )
  
} else {
  
  models <- list(
    
    # 1: NULL
    
    null = as.formula("stopover ~ 1"),
    
    # 2: DROUGHT — additive only
    
    drought = as.formula("stopover ~
                         s(SPEI_spring_1995_2020) +
                         s(SPEI_prev_spring_yr) +
                         s(SPEI_prev_fall_szn) +
                         s(SPEI_prev_winter) +
                         s(year, bs='re') +
                         s(radar_id, bs='re')"),
    
    # 3: CLIMATE — additive only
    
    climate_additive = as.formula("stopover ~
                         s(SPEI_spring_1995_2020) +
                         s(SPEI_prev_spring_yr) +
                         s(SPEI_prev_fall_szn) +
                         s(SPEI_prev_winter) +
                         s(temp_spring_PRISM_stack) +
                         s(ppt_spring_PRISM_stack) +
                         s(uwind_spring, k=20) +
                         s(vwind_spring, k=20) +
                         s(year, bs='re') +
                         s(radar_id, bs='re')"),
    
    # 4: CLIMATE — with interactions
    
    climate_interaction = as.formula("stopover ~
                         s(SPEI_spring_1995_2020) +
                         s(SPEI_prev_spring_yr) +
                         s(SPEI_prev_fall_szn) +
                         s(SPEI_prev_winter) +
                         s(temp_spring_PRISM_stack) +
                         s(ppt_spring_PRISM_stack) +
                         s(uwind_spring, k=20) +
                         s(vwind_spring, k=20) +
                         ti(SPEI_spring_1995_2020, SPEI_prev_winter, k=c(9,9)) +
                         ti(SPEI_spring_1995_2020, temp_spring_PRISM_stack, k=c(9,9)) +
                         s(year, bs='re') +
                         s(radar_id, bs='re')"),
    
    # 5: HABITAT — additive only
    
    habitat_additive = as.formula("stopover ~
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
                         s(year, bs='re') +
                         s(radar_id, bs='re')"),
    
    # 6: HABITAT — with interactions
    
    habitat_interaction = as.formula("stopover ~
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
                         ti(SPEI_spring_1995_2020, SPEI_prev_winter, k=c(9,9)) +
                         s(SPEI_spring_1995_2020, NLCD_1995_2020, bs='fs', k=4) +
                         s(SPEI_spring_1995_2020, ecoregion, bs='fs', k=4) +
                         s(SPEI_prev_winter, ecoregion, bs='fs', k=4) +
                         ti(SPEI_spring_1995_2020, alan_dist_stack, k=c(9,9)) +
                         ti(SPEI_prev_winter, dwater_temp_stack, k=c(9,9)) +
                         ti(SPEI_spring_1995_2020, dwater_temp_stack, k=c(9,9)) +
                         ti(SPEI_prev_winter, msavi_spring, k=c(9,9)) +
                         s(year, bs='re') +
                         s(radar_id, bs='re')"),
    
    # 7: DROUGHT INTERACTION — drought terms + concurrent x previous season
    
    drought_interaction = as.formula("stopover ~
                         s(SPEI_spring_1995_2020) +
                         s(SPEI_prev_spring_yr) +
                         s(SPEI_prev_fall_szn) +
                         s(SPEI_prev_winter) +
                         ti(SPEI_spring_1995_2020, SPEI_prev_winter, k=c(9,9)) +
                         s(year, bs='re') +
                         s(radar_id, bs='re')"),
    
    # 8: FULL — additive only
    
    full_additive = as.formula("stopover ~
                         s(SPEI_spring_1995_2020) +
                         s(SPEI_prev_spring_yr) +
                         s(SPEI_prev_fall_szn) +
                         s(SPEI_prev_winter) +
                         s(temp_spring_PRISM_stack) +
                         s(ppt_spring_PRISM_stack) +
                         s(uwind_spring, k=20) +
                         s(vwind_spring, k=20) +
                         s(alan_dist_stack) +
                         s(dcoast) +
                         s(dwater_perm_stack) +
                         s(dwater_temp_stack) +
                         s(forest_age_stack) +
                         s(msavi_spring) +
                         NLCD_1995_2020 +
                         s(year, bs='re') +
                         s(radar_id, bs='re')"),
    
    # 9: FULL — with all interactions
    
    full_interaction = as.formula("stopover ~
                         s(SPEI_spring_1995_2020) +
                         s(SPEI_prev_spring_yr) +
                         s(SPEI_prev_fall_szn) +
                         s(SPEI_prev_winter) +
                         s(temp_spring_PRISM_stack) +
                         s(ppt_spring_PRISM_stack) +
                         s(uwind_spring, k=20) +
                         s(vwind_spring, k=20) +
                         s(alan_dist_stack) +
                         s(dcoast) +
                         s(dwater_perm_stack) +
                         s(dwater_temp_stack) +
                         s(forest_age_stack) +
                         s(msavi_spring) +
                         NLCD_1995_2020 +
                         ti(SPEI_spring_1995_2020, SPEI_prev_winter, k=c(9,9)) +
                         s(SPEI_spring_1995_2020, NLCD_1995_2020, bs='fs', k=4) +
                         s(SPEI_spring_1995_2020, ecoregion, bs='fs', k=4) +
                         s(SPEI_prev_winter, ecoregion, bs='fs', k=4) +
                         ti(SPEI_spring_1995_2020, alan_dist_stack, k=c(9,9)) +
                         ti(uwind_spring, alan_dist_stack, k=c(9,9)) +
                         ti(ppt_spring_PRISM_stack, msavi_spring, k=c(9,9)) +
                         ti(SPEI_spring_1995_2020, temp_spring_PRISM_stack, k=c(9,9)) +
                         ti(SPEI_prev_winter, dwater_temp_stack, k=c(9,9)) +
                         ti(SPEI_spring_1995_2020, dwater_temp_stack, k=c(9,9)) +
                         ti(SPEI_prev_winter, msavi_spring, k=c(9,9)) +
                         s(year, bs='re') +
                         s(radar_id, bs='re')")
  )
}

cat("Total models defined:", length(models), "\n")
for(m in names(models)) cat(" -", m, "\n")

# PARALLEL MODEL FITTING ------------------------------------

cat("\n", paste(rep("=", 50), collapse=""), "\n")
cat("FITTING", length(models), "MODELS IN PARALLEL\n")
cat(paste(rep("=", 50), collapse=""), "\n\n")

# Determine number of workers and threads per model based on available cores

n_workers         <- min(n_cores, length(models))
threads_per_model <- max(1, floor(n_cores / length(models)))

cat("Workers:", n_workers, "\n")
cat("Threads per model:", threads_per_model, "\n")

# Set up parallel backend

cl <- makeCluster(n_workers)
registerDoParallel(cl)
clusterExport(cl, c("pooled_data", "models", "threads_per_model"),
              envir = environment())
clusterEvalQ(cl, { library(mgcv); library(dplyr) })

# Tracking

model_fit_start <- Sys.time()

# Fit each model in parallel using foreach

results <- foreach(
  model_name     = names(models),
  .packages      = c('mgcv', 'dplyr'),
  .errorhandling = 'pass',
  .verbose       = FALSE
) %dopar% {
  
  model_start <- Sys.time()
  
  tryCatch({
    
    if(model_name == "null") {
      fitted_model <- gam(
        formula = models[[model_name]],
        data    = pooled_data,
        family  = tw(link = "log"),
        method  = "REML"
      )
    } else {
      fitted_model <- bam(
        formula    = models[[model_name]],
        data       = pooled_data,
        family     = tw(link = "log"),
        method     = "fREML",
        discrete   = TRUE,
        nthreads   = threads_per_model,
        gc.level   = 1,
        control    = gam.control(trace = FALSE, epsilon = 1e-6, maxit = 200),
        samfrac    = 0.2,
        chunk.size = 10000
      )
    }
    
    fitting_time  <- as.numeric(difftime(Sys.time(), model_start, units = "mins"))
    model_summary <- summary(fitted_model)
    
    list(
      model   = fitted_model,
      stats   = data.frame(
        model         = model_name,
        aic           = AIC(fitted_model),
        r_squared     = model_summary$r.sq,
        dev_explained = model_summary$dev.expl,
        edf           = sum(fitted_model$edf),
        fitting_time  = fitting_time,
        stringsAsFactors = FALSE
      ),
      success = TRUE
    )
    
  }, error = function(e) {
    list(
      model   = NULL,
      stats   = data.frame(
        model         = model_name,
        aic           = NA,
        r_squared     = NA,
        dev_explained = NA,
        edf           = NA,
        fitting_time  = as.numeric(difftime(Sys.time(), model_start,
                                            units = "mins")),
        stringsAsFactors = FALSE
      ),
      success = FALSE,
      error   = e$message
    )
  })
}

model_fit_time <- difftime(Sys.time(), model_fit_start, units = "mins")
cat("\nAll models fitted in", round(model_fit_time, 2), "minutes\n")

stopCluster(cl)

# PROCESS RESULTS ------------------------------------

cat("\n", paste(rep("=", 50), collapse=""), "\n")
cat("PROCESSING RESULTS\n")
cat(paste(rep("=", 50), collapse=""), "\n\n")

fitted_models <- list()
model_results <- data.frame()

# Iterate through results and extract model objects and stats

for(i in seq_along(results)) {
  model_name <- names(models)[i]
  
  if(!is.list(results[[i]]) || is.null(results[[i]]$success)) {
    cat("✗", toupper(model_name), "- FAILED (foreach error)\n")
    model_results <- rbind(model_results, data.frame(
      model=model_name, aic=NA, r_squared=NA,
      dev_explained=NA, edf=NA, fitting_time=NA,
      stringsAsFactors=FALSE))
    next
  }
  
  if(results[[i]]$success) {
    fitted_models[[model_name]] <- results[[i]]$model
    model_results <- rbind(model_results, results[[i]]$stats)
    cat("✓", toupper(model_name),
        "- AIC:", round(results[[i]]$stats$aic, 2),
        "R²:", round(results[[i]]$stats$r_squared, 3),
        "Dev:", round(results[[i]]$stats$dev_explained, 3),
        "Time:", round(results[[i]]$stats$fitting_time, 2), "min\n")
  } else {
    cat("✗", toupper(model_name), "- FAILED:", results[[i]]$error, "\n")
    model_results <- rbind(model_results, results[[i]]$stats)
  }
}

if(nrow(model_results) == 0) stop("All models failed — check log for errors")

# MODEL COMPARISON ------------------------------------

# Rank models by AIC, calculate delta AIC and AIC weights for comparison.

model_results           <- model_results[order(model_results$aic, na.last = TRUE), ]
model_results$delta_aic <- model_results$aic - min(model_results$aic, na.rm = TRUE)
model_results$aic_weight <- exp(-0.5 * model_results$delta_aic) /
  sum(exp(-0.5 * model_results$delta_aic), na.rm = TRUE)

cat("\nMODEL COMPARISON TABLE:\n")
print(model_results[, c("model", "aic", "delta_aic", "aic_weight",
                        "r_squared", "dev_explained", "fitting_time")])

best_model_name <- model_results$model[1]
cat("\nBEST MODEL:", toupper(best_model_name), "\n")
cat("AIC:", round(model_results$aic[1], 2), "\n")
cat("Deviance explained:", round(model_results$dev_explained[1], 3), "\n")

# SAVE RESULTS ------------------------------------

cat("\nSaving results...\n")

saveRDS(fitted_models,
        file.path(output_dir, paste0(season, "_all_models.rds")))

write.csv(model_results,
          file.path(output_dir, paste0(season, "_model_comparison.csv")),
          row.names = FALSE)

if(best_model_name %in% names(fitted_models)) {
  saveRDS(fitted_models[[best_model_name]],
          file.path(output_dir, paste0(season, "_best_model.rds")))
}

write.csv(pooled_data,
          file.path(output_dir, paste0(season, "_gam_data.csv")),
          row.names = FALSE)

total_time <- difftime(Sys.time(), total_start, units = "mins")

cat(paste(rep("=", 60), collapse=""), "\n")
cat("COMPLETE\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("Total analysis time:", round(total_time, 2), "minutes\n")
cat("Results saved to:", output_dir, "\n")