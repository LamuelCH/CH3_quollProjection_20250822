#!/usr/bin/env Rscript
# ============================================================================
# Multi-species Spatial Prediction under Climate Change Scenarios
# With SLURM Array Job Support
# ============================================================================
# Description: This script performs spatial predictions for multiple species
#              using a fitted sfMsPGOcc model under various climate change
#              scenarios, with support for parallel processing via job arrays
# 
# Author: Lamuel C.H. Chung
# Date Created: 2025-01-01
# Last Modified: 2025-01-01
# ============================================================================

# Load required libraries
cat("\n========================================\n")
cat("LOADING REQUIRED LIBRARIES\n")
cat("========================================\n")
library(terra)
library(dplyr)
library(spOccupancy)

# ============================================================================
# JOB ARRAY CONFIGURATION
# ============================================================================
cat("\n========================================\n")
cat("JOB ARRAY CONFIGURATION\n")
cat("========================================\n")

# Parse command line arguments for SLURM array jobs
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
  job_index <- 0  # Default for testing
  cat("No job index provided, using default value 0 for local testing\n")
} else {
  job_index <- as.integer(args[1])
  if(is.na(job_index)) {
    stop("Invalid job index: must be an integer")
  }
  cat("Job index from SLURM:", job_index, "\n")
}

# Define all scenario combinations
models <- c("gfdl-esm4", "mpi-esm1-2-hr", "ukesm1-0-ll", "mri-esm2-0", "ipsl-cm6a-lr")
ssps <- c("ssp585", "ssp126", "ssp370")
time_periods <- c("2071-2100", "2011-2040", "2041-2070")

# Create all combinations
scenarios <- expand.grid(
  model = models,
  ssp = ssps,
  period = time_periods,
  stringsAsFactors = FALSE
)

# Add baseline scenario as index 0
if(job_index == 0) {
  # Baseline scenario (contemporary climate)
  current_scenario <- list(
    model = "baseline",
    ssp = "current",
    period = "1981-2010"
  )
  is_baseline <- TRUE
  cat("\n=== BASELINE SCENARIO ===\n")
  cat("Running prediction for contemporary climate (1981-2010)\n")
} else {
  # Validate job index
  if(job_index > nrow(scenarios)) {
    stop(sprintf("Job index %d out of range (max: %d)", job_index, nrow(scenarios)))
  }
  
  # Get current scenario
  current_scenario <- scenarios[job_index, ]
  is_baseline <- FALSE
  
  cat("\n=== CLIMATE CHANGE SCENARIO ===\n")
  cat("Scenario", job_index, "of", nrow(scenarios), "\n")
  cat("  Model:", current_scenario$model, "\n")
  cat("  SSP:", current_scenario$ssp, "\n")
  cat("  Period:", current_scenario$period, "\n")
}

# ============================================================================
# CONFIGURATION
# ============================================================================
cat("\n========================================\n")
cat("CONFIGURATION\n")
cat("========================================\n")

# Set parameters
CHUNK_SIZE <- 1000
N_THREADS <- 8
SAVE_INTERVAL <- 5  # Save intermediate results every N chunks

# Create scenario-specific output directory
if(is_baseline) {
  scenario_suffix <- "baseline_1981-2010"
} else {
  scenario_suffix <- paste(current_scenario$period, current_scenario$model, 
                           current_scenario$ssp, sep = "_")
}

current_date <- format(Sys.Date(), "%Y%m%d")
base_output_dir <- paste0("output/mainland/sfMsPGOcc/predictions/climate_projections/",
                          scenario_suffix, "/pred_", current_date)
chunks_dir <- file.path(base_output_dir, "chunks")

cat("Output directory:", base_output_dir, "\n")
cat("Scenario identifier:", scenario_suffix, "\n")
cat("Chunk size:", CHUNK_SIZE, "locations per chunk\n")
cat("Number of threads:", N_THREADS, "\n")

# Create output directories
dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(chunks_dir, recursive = TRUE, showWarnings = FALSE)
cat("✓ Output directories created\n")

# ============================================================================
# LOAD MODEL AND DATA
# ============================================================================
cat("\n========================================\n")
cat("LOADING MODEL AND DATA\n")
cat("========================================\n")

# Load the fitted model
model_file <- "models/mainland/sfMsPGOcc/model_20250718_sfMsPGOcc_2010-2024_nthin500_nburn5e+05_nchain1_nsample6000_nfactors3_neighbors10.RData"
cat("Loading model from:", model_file, "\n")
load(model_file)
cat("✓ Model loaded successfully\n")

# Load the original data (for standardization parameters)
data_file <- "input/mainland/data_top12_2010-2024_5r.RData"
cat("Loading original data from:", data_file, "\n")
load(data_file)
cat("✓ Data loaded successfully\n")

# Determine number of species
if (!is.null(out.sfMsPGOcc$y)) {
  n_species <- dim(out.sfMsPGOcc$y)[1]
  species_names <- rownames(out.sfMsPGOcc$y)
  cat("\nModel contains", n_species, "species:\n")
  for (i in 1:n_species) {
    cat("  Species", i, ":", species_names[i], "\n")
  }
} else {
  n_species <- 12
  species_names <- paste0("Species_", 1:n_species)
  cat("\n⚠ Could not determine species names from model\n")
  cat("  Assuming", n_species, "species\n")
}

# ============================================================================
# LOAD AND PREPARE ENVIRONMENTAL DATA
# ============================================================================
cat("\n========================================\n")
cat("LOADING ENVIRONMENTAL DATA\n")
cat("========================================\n")

# Load bioclimatic variables based on scenario
if(is_baseline) {
  # Load baseline CHELSA data
  cat("Loading baseline CHELSA bioclimatic variables (1981-2010)...\n")
  bio <- rast("input/mainland/1981-2010/mainland_CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
} else {
  # Construct path to climate projection file
  bio_file <- file.path("input/mainland", 
                        current_scenario$period,
                        current_scenario$ssp,
                        paste0("mainland_CHELSA_bio_", 
                               current_scenario$period, "_",
                               current_scenario$model, "_",
                               current_scenario$ssp, "_V.2.1_EPSG3577.tif"))
  
  cat("Loading climate projection from:", bio_file, "\n")
  
  # Check if file exists
  if(!file.exists(bio_file)) {
    stop(paste("Climate projection file not found:", bio_file))
  }
  
  bio <- rast(bio_file)
}

# Select relevant bioclimatic variables
bio <- bio[[c("bio5", "bio6", "bio15")]]
cat("  Selected variables: bio5, bio6, bio15\n")

# Load static environmental layers (terrain and vegetation)
cat("Loading terrain variables (static)...\n")
terrain <- rast("input/mainland/mainland_env_terrain_EPSG3577.tif")
terrain <- terrain[[c("dem", "tpi")]]
cat("  Selected variables: dem, tpi\n")

cat("Loading foliage projective cover (static)...\n")
fpc <- rast("input/mainland/mainland_env_fpc_EPSG3577.tif")
cat("  Variable: fpc\n")

# Stack all environmental layers
env_stack <- c(bio, terrain, fpc)
cat("\nEnvironmental stack created with", nlyr(env_stack), "layers\n")

# Print raster information
cat("\nRaster properties:\n")
cat("  CRS:", crs(env_stack, describe = TRUE)$name, "\n")
cat("  Extent:", ext(env_stack)[1], "-", ext(env_stack)[2], "(x),",
    ext(env_stack)[3], "-", ext(env_stack)[4], "(y)\n")
cat("  Resolution:", res(env_stack)[1], "x", res(env_stack)[2], "meters\n")
cat("  Dimensions:", nrow(env_stack), "rows x", ncol(env_stack), "cols\n")

# Print climate scenario summary
if(!is_baseline) {
  cat("\nClimate scenario details:\n")
  cat("  Climate model:", current_scenario$model, "\n")
  cat("  SSP scenario:", current_scenario$ssp, "\n")
  cat("  Time period:", current_scenario$period, "\n")
}

# ============================================================================
# PREPARE PREDICTION DATA
# ============================================================================
cat("\n========================================\n")
cat("PREPARING PREDICTION DATA\n")
cat("========================================\n")

# Convert raster to data frame
cat("Converting raster stack to data frame...\n")
start_time <- Sys.time()
env_stack_df <- as.data.frame(env_stack, xy = TRUE, na.rm = TRUE)
cat("  Conversion completed in", round(difftime(Sys.time(), start_time, units = "secs"), 2), "seconds\n")
cat("  Total non-NA cells:", format(nrow(env_stack_df), big.mark = ","), "\n")

# Store original coordinates in meters for later raster creation
env_stack_df$x_meters <- env_stack_df$x
env_stack_df$y_meters <- env_stack_df$y

# Convert coordinates to kilometers for model prediction
cat("\nConverting coordinates to kilometers for model...\n")
env_stack_df$x <- env_stack_df$x / 1000
env_stack_df$y <- env_stack_df$y / 1000
cat("  Coordinates converted to km\n")

# Extract standardization parameters from original data
cat("\nExtracting standardization parameters from training data...\n")
orig_means <- sapply(1:ncol(data$occ.covs), function(i) mean(data$occ.covs[, i]))
orig_sds <- sapply(1:ncol(data$occ.covs), function(i) sd(data$occ.covs[, i]))

# Variable names
env_vars <- c("bio5", "bio6", "bio15", "dem", "tpi", "fpc")
std_vars <- paste0(env_vars, "_std")

# Print standardization parameters
cat("\nStandardization parameters (from training data):\n")
cat("Variable | Mean      | SD\n")
cat("---------|-----------|----------\n")
for (i in 1:length(env_vars)) {
  cat(sprintf("%-8s | %9.3f | %9.3f\n", env_vars[i], orig_means[i], orig_sds[i]))
}

# Standardize variables using training data parameters
cat("\nStandardizing environmental variables...\n")
for (i in 1:length(env_vars)) {
  env_stack_df[[std_vars[i]]] <- (env_stack_df[[env_vars[i]]] - orig_means[i]) / orig_sds[i]
  
  # Check for issues
  n_nan <- sum(is.nan(env_stack_df[[std_vars[i]]]))
  n_inf <- sum(is.infinite(env_stack_df[[std_vars[i]]]))
  
  if (n_nan > 0 || n_inf > 0) {
    cat("  ⚠", env_vars[i], ": Found", n_nan, "NaN and", n_inf, "Inf values\n")
    env_stack_df[[std_vars[i]]][is.nan(env_stack_df[[std_vars[i]]])] <- NA
    env_stack_df[[std_vars[i]]][is.infinite(env_stack_df[[std_vars[i]]])] <- NA
  } else {
    cat("  ✓", env_vars[i], "standardized successfully\n")
  }
}

# ============================================================================
# CREATE DESIGN MATRIX
# ============================================================================
cat("\n========================================\n")
cat("CREATING DESIGN MATRIX\n")
cat("========================================\n")

cat("Building design matrix with quadratic terms...\n")
X.0 <- cbind(
  1,  # Intercept
  env_stack_df$bio5_std, env_stack_df$bio5_std^2,
  env_stack_df$bio6_std, env_stack_df$bio6_std^2,
  env_stack_df$bio15_std,
  env_stack_df$dem_std, env_stack_df$dem_std^2,
  env_stack_df$tpi_std, env_stack_df$tpi_std^2,
  env_stack_df$fpc_std, env_stack_df$fpc_std^2
)

# Set column names for clarity
colnames(X.0) <- c("intercept", 
                   "bio5", "bio5_sq",
                   "bio6", "bio6_sq",
                   "bio15",
                   "dem", "dem_sq",
                   "tpi", "tpi_sq",
                   "fpc", "fpc_sq")

cat("  Design matrix created with", ncol(X.0), "predictors\n")

# Prepare coordinates in km for prediction
coords.0 <- as.matrix(env_stack_df[, c('x', 'y')])

# Keep track of original meter coordinates for raster creation
coords_meters <- as.matrix(env_stack_df[, c('x_meters', 'y_meters')])

# Check for and remove NA values
cat("\nChecking for NA values...\n")
na_rows <- rowSums(is.na(X.0)) > 0
n_na <- sum(na_rows)

if (n_na > 0) {
  cat("  ⚠ Found", format(n_na, big.mark = ","), "rows with NA values\n")
  X.0 <- X.0[!na_rows, ]
  coords.0 <- coords.0[!na_rows, ]
  coords_meters <- coords_meters[!na_rows, ]
  cat("  After removal:", format(nrow(X.0), big.mark = ","), "locations remain\n")
} else {
  cat("  ✓ No NA values found\n")
}

# ============================================================================
# CHUNK PROCESSING SETUP
# ============================================================================
cat("\n========================================\n")
cat("CHUNK PROCESSING SETUP\n")
cat("========================================\n")

n_locations <- nrow(X.0)
n_chunks <- ceiling(n_locations / CHUNK_SIZE)

cat("Total locations to predict:", format(n_locations, big.mark = ","), "\n")
cat("Chunk size:", format(CHUNK_SIZE, big.mark = ","), "\n")
cat("Number of chunks:", n_chunks, "\n")
cat("Estimated processing time: ~", round(n_chunks * 2, 1), "minutes (rough estimate)\n")

# Initialize results data frame with meter coordinates
results_summary <- data.frame(
  x = numeric(0),  # This will store coordinates in meters
  y = numeric(0),  # This will store coordinates in meters
  species = integer(0),
  species_name = character(0),
  psi_mean = numeric(0),
  psi_sd = numeric(0),
  z_prob = numeric(0),
  w_mean = numeric(0),
  w_sd = numeric(0)
)

# ============================================================================
# MAIN PREDICTION LOOP
# ============================================================================
cat("\n========================================\n")
cat("STARTING SPATIAL PREDICTIONS\n")
cat("========================================\n")

overall_start_time <- Sys.time()

for (i in 1:n_chunks) {
  chunk_start_time <- Sys.time()
  
  cat("\n----------------------------------------\n")
  cat("CHUNK", i, "OF", n_chunks, "\n")
  cat("----------------------------------------\n")
  
  # Define indices for this chunk
  start_idx <- (i - 1) * CHUNK_SIZE + 1
  end_idx <- min(i * CHUNK_SIZE, n_locations)
  chunk_size_actual <- end_idx - start_idx + 1
  
  cat("Processing locations", format(start_idx, big.mark = ","), 
      "to", format(end_idx, big.mark = ","),
      "(", chunk_size_actual, "locations)\n")
  
  # Extract chunk data
  X.0.chunk <- X.0[start_idx:end_idx, , drop = FALSE]
  coords.0.chunk <- coords.0[start_idx:end_idx, , drop = FALSE]  # km coordinates for prediction
  coords_meters.chunk <- coords_meters[start_idx:end_idx, , drop = FALSE]  # meter coordinates for output
  
  # Run prediction with error handling
  tryCatch({
    cat("Running spatial predictions...\n")
    
    pred.chunk <- predict(out.sfMsPGOcc, 
                          X.0.chunk,  
                          coords.0.chunk,  # Use km coordinates for prediction
                          n.omp.threads = N_THREADS,
                          verbose = TRUE,
                          n.report = 1000,
                          type = 'occupancy')
    
    cat("✓ Predictions completed\n")
    
    # Process results
    cat("Processing prediction results...\n")
    chunk_summary <- data.frame()
    
    # Get dimensions
    n_loc <- dim(pred.chunk$psi.0.samples)[3]
    w_dims <- dim(pred.chunk$w.0.samples)
    has_species_w <- length(w_dims) >= 3 && w_dims[2] > 1
    
    if (has_species_w) {
      cat("  - Spatial random effects: Species-specific\n")
    } else {
      cat("  - Spatial random effects: Shared across species\n")
    }
    
    # Process each location and species
    for (loc in 1:n_loc) {
      # Extract spatial random effect for this location
      if (has_species_w) {
        w_means <- apply(pred.chunk$w.0.samples[, , loc], 2, mean)
        w_sds <- apply(pred.chunk$w.0.samples[, , loc], 2, sd)
      } else {
        w_mean <- mean(pred.chunk$w.0.samples[, 1, loc])
        w_sd <- sd(pred.chunk$w.0.samples[, 1, loc])
      }
      
      for (sp in 1:n_species) {
        # Calculate summary statistics
        psi_mean <- mean(pred.chunk$psi.0.samples[, sp, loc])
        psi_sd <- sd(pred.chunk$psi.0.samples[, sp, loc])
        z_prob <- mean(pred.chunk$z.0.samples[, sp, loc])
        
        # Get spatial random effect
        if (has_species_w) {
          w_mean_sp <- w_means[sp]
          w_sd_sp <- w_sds[sp]
        } else {
          w_mean_sp <- w_mean
          w_sd_sp <- w_sd
        }
        
        # Create row with METER coordinates for output
        row_data <- data.frame(
          x = coords_meters.chunk[loc, 1],  # Store in meters
          y = coords_meters.chunk[loc, 2],  # Store in meters
          species = sp,
          species_name = species_names[sp],
          psi_mean = psi_mean,
          psi_sd = psi_sd,
          z_prob = z_prob,
          w_mean = w_mean_sp,
          w_sd = w_sd_sp
        )
        
        chunk_summary <- rbind(chunk_summary, row_data)
      }
    }
    
    # Append to overall results
    results_summary <- rbind(results_summary, chunk_summary)
    
    # Save chunk results
    chunk_file <- file.path(chunks_dir, paste0("chunk_", sprintf("%03d", i), "_of_", n_chunks, ".RData"))
    save(chunk_summary, file = chunk_file)
    cat("✓ Chunk results saved to:", basename(chunk_file), "\n")
    
    # Calculate chunk statistics
    chunk_time <- difftime(Sys.time(), chunk_start_time, units = "secs")
    cat("\nChunk statistics:\n")
    cat("  - Processing time:", round(chunk_time, 2), "seconds\n")
    cat("  - Locations processed:", chunk_size_actual, "\n")
    cat("  - Total predictions:", nrow(chunk_summary), "\n")
    cat("  - Mean occurrence probability:", round(mean(chunk_summary$z_prob), 3), "\n")
    
    # Clean up memory
    rm(pred.chunk, chunk_summary)
    gc(verbose = FALSE)
    
  }, error = function(e) {
    cat("\n⚠ ERROR in chunk", i, ":", conditionMessage(e), "\n")
    
    # Implement fallback for memory errors
    if (grepl("cannot allocate", conditionMessage(e))) {
      cat("⚠ Memory allocation error detected. Consider reducing chunk size.\n")
    }
  })
  
  # Save intermediate results periodically
  if (i %% SAVE_INTERVAL == 0 || i == n_chunks) {
    intermediate_file <- file.path(base_output_dir, "intermediate_prediction_results.RData")
    save(results_summary, file = intermediate_file)
    cat("\n✓ Intermediate results saved (", format(nrow(results_summary), big.mark = ","), 
        "total predictions)\n")
  }
  
  # Progress summary
  progress_pct <- round(i / n_chunks * 100, 1)
  elapsed_time <- difftime(Sys.time(), overall_start_time, units = "mins")
  estimated_remaining <- elapsed_time / i * (n_chunks - i)
  
  cat("\nOverall progress:", progress_pct, "%\n")
  cat("Elapsed time:", round(elapsed_time, 1), "minutes\n")
  if (i < n_chunks) {
    cat("Estimated remaining:", round(estimated_remaining, 1), "minutes\n")
  }
}

# ============================================================================
# SAVE FINAL RESULTS
# ============================================================================
cat("\n========================================\n")
cat("SAVING FINAL RESULTS\n")
cat("========================================\n")

# Save combined results with scenario information
final_results_file <- file.path(base_output_dir, paste0("combined_prediction_results_", 
                                                        scenario_suffix, ".RData"))
save(results_summary, current_scenario, file = final_results_file)
cat("✓ Final results saved to:", final_results_file, "\n")
cat("  Total predictions:", format(nrow(results_summary), big.mark = ","), "\n")

# ============================================================================
# CREATE SPECIES RASTERS
# ============================================================================
cat("\n========================================\n")
cat("CREATING SPECIES RASTERS\n")
cat("========================================\n")

if (nrow(results_summary) > 0) {
  species_ids <- unique(results_summary$species)
  
  for (sp in species_ids) {
    cat("\nProcessing rasters for", species_names[sp], "(Species", sp, ")...\n")
    
    # Filter data for this species
    sp_data <- results_summary[results_summary$species == sp, ]
    
    # Create species output directory with scenario suffix
    species_dir <- file.path(base_output_dir, 
                             paste0("species_", sprintf("%02d", sp), "_", 
                                    species_names[sp], "_", scenario_suffix))
    dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Create rasters using METER coordinates
    cat("  Creating rasters (coordinates in meters)...\n")
    
    # Occurrence probability
    psi_mean_raster <- rast(sp_data[, c("x", "y", "psi_mean")], 
                            type = "xyz", crs = "EPSG:3577")
    writeRaster(psi_mean_raster, 
                file.path(species_dir, paste0("psi_mean_", scenario_suffix, ".tif")), 
                overwrite = TRUE)
    cat("    ✓ psi_mean.tif (mean occurrence probability)\n")
    
    # Uncertainty
    psi_sd_raster <- rast(sp_data[, c("x", "y", "psi_sd")], 
                          type = "xyz", crs = "EPSG:3577")
    writeRaster(psi_sd_raster, 
                file.path(species_dir, paste0("psi_sd_", scenario_suffix, ".tif")), 
                overwrite = TRUE)
    cat("    ✓ psi_sd.tif (occurrence probability SD)\n")
    
    # Posterior occurrence probability
    z_prob_raster <- rast(sp_data[, c("x", "y", "z_prob")], 
                          type = "xyz", crs = "EPSG:3577")
    writeRaster(z_prob_raster, 
                file.path(species_dir, paste0("z_prob_", scenario_suffix, ".tif")), 
                overwrite = TRUE)
    cat("    ✓ z_prob.tif (posterior occurrence probability)\n")
    
    # Spatial random effects
    w_mean_raster <- rast(sp_data[, c("x", "y", "w_mean")], 
                          type = "xyz", crs = "EPSG:3577")
    writeRaster(w_mean_raster, 
                file.path(species_dir, paste0("w_mean_", scenario_suffix, ".tif")), 
                overwrite = TRUE)
    cat("    ✓ w_mean.tif (mean spatial random effect)\n")
    
    w_sd_raster <- rast(sp_data[, c("x", "y", "w_sd")], 
                        type = "xyz", crs = "EPSG:3577")
    writeRaster(w_sd_raster, 
                file.path(species_dir, paste0("w_sd_", scenario_suffix, ".tif")), 
                overwrite = TRUE)
    cat("    ✓ w_sd.tif (spatial random effect SD)\n")
    
    # Summary statistics
    cat("  Summary statistics:\n")
    cat("    - Mean occurrence probability:", round(mean(sp_data$z_prob), 3), "\n")
    cat("    - SD occurrence probability:", round(sd(sp_data$z_prob), 3), "\n")
    cat("    - Min/Max occurrence:", round(min(sp_data$z_prob), 3), "/", 
        round(max(sp_data$z_prob), 3), "\n")
  }
}

# ============================================================================
# SAVE SCENARIO METADATA
# ============================================================================
cat("\n========================================\n")
cat("SAVING SCENARIO METADATA\n")
cat("========================================\n")

metadata <- list(
  job_index = job_index,
  scenario = current_scenario,
  is_baseline = is_baseline,
  scenario_suffix = scenario_suffix,
  processing_date = Sys.Date(),
  processing_time = difftime(Sys.time(), overall_start_time, units = "mins"),
  n_species = n_species,
  n_locations = n_locations,
  n_chunks = n_chunks,
  chunk_size = CHUNK_SIZE,
  model_file = model_file,
  output_dir = base_output_dir
)

metadata_file <- file.path(base_output_dir, paste0("metadata_", scenario_suffix, ".RData"))
save(metadata, file = metadata_file)
cat("✓ Metadata saved to:", metadata_file, "\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("PREDICTION COMPLETE\n")
cat("========================================\n")

total_time <- difftime(Sys.time(), overall_start_time, units = "mins")
cat("Total processing time:", round(total_time, 2), "minutes\n")
cat("Output directory:", base_output_dir, "\n")
cat("Scenario:", scenario_suffix, "\n")

if (!is_baseline) {
  cat("\nClimate projection details:\n")
  cat("  - Model:", current_scenario$model, "\n")
  cat("  - SSP:", current_scenario$ssp, "\n")
  cat("  - Period:", current_scenario$period, "\n")
}

cat("\nFiles created:\n")
cat("  - Combined results:", paste0("combined_prediction_results_", scenario_suffix, ".RData\n"))
cat("  - Chunk results:", n_chunks, "files in chunks/\n")
cat("  - Species rasters:", n_species, "folders with 5 rasters each\n")
cat("  - Metadata:", paste0("metadata_", scenario_suffix, ".RData\n"))

cat("\n✓ Multi-species spatial prediction completed successfully!\n")
cat("========================================\n\n")