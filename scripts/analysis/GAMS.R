##############################
#
#   GAMS - Fall & Spring 
#   Ian Becker
#   December 2025
#
##############################

# This script is used to run GAMs for both fall and spring seasons

library(mgcv)
library(parallel)
library(doParallel)
library(foreach)

setwd("PATH HERE")

#################################
#   PREP DATA
#################################

cat("\nPREPARING DATA FOR GAMs\n")

# Select season (change between fall and spring as needed)

season <- "fall"  

# Set paths

input_data <- paste0(season, "_gam_data.csv")  # Pre-processed dataset
output_dir <- paste0("gam_results/", season, "_model_comparison")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Parallel processing

n_cores <- XXXX

# Read in data

cat("LOADING DATA\n")
pooled_data <- read.csv(input_data, stringsAsFactors = FALSE)

# Data structure checks

cat(toupper(season), "DATA:")
cat("OBSERVATIONS:", nrow(pooled_data), "\n")
cat("VARIABLES:", ncol(pooled_data), "\n")

# Convert categorical variables to factors for model input

categorical_vars <- c("ecoregion", "NLCD_1995_2020")
for(var in categorical_vars) {
  if(var %in% names(pooled_data)) {
    pooled_data[[var]] <- as.factor(pooled_data[[var]])
    pooled_data[[var]] <- droplevels(pooled_data[[var]])
    cat("CONVERTED", var, "TO", nlevels(pooled_data[[var]]), "FACTOR LEVELS\n")
  }
}

cat("DATA PREP COMPLETE\n")

#################################
#   SETUP MODELS
#################################

cat("\nSETTING UP MODELS FOR", toupper(season), "\n")

# If else statement to define models based on season

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

cat("DEFINED", length(models), "FOR", toupper(season),"\n")

#################################
# SETUP PARALLEL PROCESSING
#################################

cat("\nSETTING UP PARALLEL PROCESSING\n")

# Create cluster - use one worker per model (up to n_cores)

cl <- makeCluster(min(n_cores, length(models)))

# Calculate threads per model

threads_per_model <- max(1, floor(n_cores / length(models)))
cat("USING", length(cl), "WORKERS WITH", threads_per_model, "THREADS EACH\n")

# Register the cluster as the parallel backend for foreach

registerDoParallel(cl)

# Verify registration

cat("PARALLEL BACKED REGISTERED", getDoParName(), "\n")
cat("# OF WORKERS", getDoParWorkers(), "\n")

# Export necessary objects to ALL workers in the cluster

clusterExport(cl, c("pooled_data", "models", "threads_per_model"), envir = environment())

# Load required packages on ALL workers

clusterEvalQ(cl, {
  library(mgcv)
})

cat("PARALLEL WORKERS INTIALIZED\n")

#################################
#  FITTING MODELS
#################################

cat("\nFITTING", length(models), "MODELS IN PARALLEL\n")

# Track timing for each model

model_fit_start <- Sys.time()

# Fitting models in parallel 

results <- foreach(
  model_name = names(models),
  .packages = 'mgcv',
  .errorhandling = 'pass',  # Continues even if one model fails
  .verbose = FALSE
) %dopar% {  # %dopar% uses the registered parallel backend (cl)
  
  model_start <- Sys.time()
  
  tryCatch({
    
    # Fitting null model separately since it requires simpler arguments
    
    if(model_name == "null") {
      fitted_model <- gam(
        formula = models[[model_name]],
        data = pooled_data,
        family = nb(),
        method = "REML"
      )
      
    } else {
      
      # Complex models - optimized for larger models
      
      fitted_model <- bam(
        formula = models[[model_name]],
        data = pooled_data,
        family = nb(),
        method = "fREML",        
        discrete = TRUE,        
        nthreads = threads_per_model,  # Auto-calculated based on n_cores
        gc.level = 1,           # Moderate cleaning for optimization
        control = gam.control(
          trace = FALSE,
          epsilon = 1e-6,
          maxit = 200
        ),
        samfrac = 0.1,          # Use 10% subsample for smoothness selection
        chunk.size = 10000      # Process data in 10k chunks
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
    
    # Return error information if neccessary
    
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
cat("MODELS FIT IN", round(total_fit_time, 2), "MINUTES\n")

#################################
# PROCESS MODEL RESULTS
#################################

cat("\nPROCESSING RESULTS\n")

# Creating data frame with all model results

fitted_models <- list()
model_results <- data.frame()

for(i in seq_along(results)) {
  model_name <- names(models)[i]
  
  if(results[[i]]$success) {
    fitted_models[[model_name]] <- results[[i]]$model
    model_results <- rbind(model_results, results[[i]]$stats)
  } else {
    model_results <- rbind(model_results, results[[i]]$stats)
  }
}

# Stop the cluster when done

stopCluster(cl)
cat("ANALYSIS COMPLETE; PARALLEL STRUCTURE STOPPED\n")

#################################
# MODEL COMPARISON 
#################################

cat("\nMODEL COMPARISON RESULTS\n")

# Sort by AIC (best model first)

model_results <- model_results[order(model_results$aic, na.last = TRUE), ]

# Calculate delta AIC

model_results$delta_aic <- model_results$aic - min(model_results$aic, na.rm = TRUE)

# Calculate AIC weights

model_results$aic_weight <- exp(-0.5 * model_results$delta_aic) / sum(exp(-0.5 * model_results$delta_aic), na.rm = TRUE)

# Print results table

cat("MODEL COMPARISON TABLE (sorted by AIC):\n")
print(model_results[, c("model", "aic", "delta_aic", "aic_weight", "r_squared", "dev_explained", "edf", "fitting_time")])

#################################
# SAVE RESULTS
#################################

cat("\nSAVING RESULTS\n")

# Save all fitted models

models_file <- file.path(output_dir, paste0(season, "_all_models.rds"))
saveRDS(fitted_models, models_file)
cat("ALL MODELS SAVED", models_file)

# Save comparison results

comparison_file <- file.path(output_dir, paste0(season, "_model_comparison.csv"))
write.csv(model_results, comparison_file, row.names = FALSE)
cat("MODEL COMPARISON SAVED", comparison_file)

# Save best model separately

best_model_file <- file.path(output_dir, paste0(season, "_best_model.rds"))
saveRDS(fitted_models[[best_model]], best_model_file)
cat("TOP MODEL SAVED", best_model_file, "\n")

cat("\nGAMS COMPLETE! RESULTS SAVED TO:", output_dir, "\n")