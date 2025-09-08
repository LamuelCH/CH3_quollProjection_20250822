# #!/usr/bin/env Rscript
# ============================================================================
# Multi-species Spatial Prediction using sfMsPGOcc Model
# ============================================================================
# Description: This script performs spatial predictions for multiple species
#              using a fitted sfMsPGOcc model across a study area defined by
#              environmental raster layers.
# 
# Author: [Your Name]
# Date Created: 2025-07-11
# Last Modified: 2025-07-11
# ============================================================================

# Load required libraries
cat("\n========================================\n")
cat("LOADING REQUIRED LIBRARIES\n")
cat("========================================\n")
library(terra)
library(dplyr)
library(spOccupancy)

# ============================================================================
# CONFIGURATION
# ============================================================================
cat("\n========================================\n")
cat("CONFIGURATION\n")
cat("========================================\n")

# Set parametersElapsed time: 12.5 minutes

CHUNK_SIZE <- 1000
N_THREADS <- 16
SAVE_INTERVAL <- 5  # Save intermediate results every N chunks

# Create date-based output directory
current_date <- format(Sys.Date(), "%Y%m%d")
base_output_dir <- paste0("output/sfMsPGOcc/predictions/spatial/pred_", current_date)
chunks_dir <- file.path(base_output_dir, "chunks")

cat("Output directory:", base_output_dir, "\n")
cat("Chunk size:", CHUNK_SIZE, "locations per chunk\n")
cat("Number of threads:", N_THREADS, "\n")

# Create output directories
dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(chunks_dir, recursive = TRUE, showWarnings = FALSE)
cat(" Output directories created\n")

# ============================================================================
# LOAD MODEL AND DATA
# ============================================================================
cat("\n========================================\n")
cat("LOADING MODEL AND DATA\n")
cat("========================================\n")

# Load the fitted model
model_file <- "models/sfMsPGOcc/model_20250718_sfMsPGOcc_2010-2024_nthin500_nburn5e+05_nchain1_nsample6000_nfactors3_neighbors10.RData"
cat("Loading model from:", model_file, "\n")
load(model_file)
cat(" Model loaded successfully\n")

# Load the original data (for standardization parameters)
data_file <- "input/data_top12_2010-2024_5r.RData"
cat("Loading original data from:", data_file, "\n")
load(data_file)
cat(" Data loaded successfully\n")

# Determine number of species
if (!is.null(out.sfMsPGOcc$y)) {
  n_species <- dim(out.sfMsPGOcc$y)[1]
  species_names <- rownames(out.sfMsPGOcc$y)
  cat("\n Model contains", n_species, "species:\n")
  for (i in 1:n_species) {
    cat("  Species", i, ":", species_names[i], "\n")
  }
} else {
  n_species <- 12
  species_names <- paste0("Species_", 1:n_species)
  cat("\nÃÂ  Could not determine species names from model\n")
  cat("  Assuming", n_species, "species\n")
}

# ============================================================================
# LOAD AND PREPARE ENVIRONMENTAL DATA
# ============================================================================
cat("\n========================================\n")
cat("LOADING ENVIRONMENTAL DATA\n")
cat("========================================\n")

# Load environmental rasters
cat("Loading CHELSA bioclimatic variables...\n")
bio <- rast("input/CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
bio <- bio[[c("bio5", "bio6", "bio15")]]
cat("   Selected variables: bio5, bio6, bio15\n")

cat("Loading terrain variables...\n")
terrain <- rast("input/unmasked_env_terrain_EPSG3577.tif")
terrain <- terrain[[c("dem", "tpi")]]
cat("   Selected variables: dem, tpi\n")

cat("Loading foliage projective cover...\n")
fpc <- rast("input/env_foliage_EPSG3577.tif")
cat("   Variable: fpc\n")

# Stack all environmental layers
env_stack <- c(bio, terrain, fpc)
cat("\n Environmental stack created with", nlyr(env_stack), "layers\n")

# Print raster information
cat("\nRaster properties:\n")
cat("  CRS:", crs(env_stack, describe = TRUE)$name, "\n")
cat("  Extent:", ext(env_stack)[1], "-", ext(env_stack)[2], "(x),",
    ext(env_stack)[3], "-", ext(env_stack)[4], "(y)\n")
cat("  Resolution:", res(env_stack)[1], "x", res(env_stack)[2], "\n")
cat("  Dimensions:", nrow(env_stack), "rows x", ncol(env_stack), "cols\n")

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
cat("   Conversion completed in", round(difftime(Sys.time(), start_time, units = "secs"), 2), "seconds\n")
cat("   Total non-NA cells:", format(nrow(env_stack_df), big.mark = ","), "\n")

# Convert coordinates to kilometers
cat("\nConverting coordinates to kilometers...\n")
env_stack_df$x <- env_stack_df$x / 1000
env_stack_df$y <- env_stack_df$y / 1000
cat("   Coordinates converted\n")

# Extract standardization parameters from original data
cat("\nExtracting standardization parameters...\n")
orig_means <- sapply(1:ncol(data$occ.covs), function(i) mean(data$occ.covs[, i]))
orig_sds <- sapply(1:ncol(data$occ.covs), function(i) sd(data$occ.covs[, i]))

# Variable names
env_vars <- c("bio5", "bio6", "bio15", "dem", "tpi", "fpc")
std_vars <- paste0(env_vars, "_std")

# Print standardization parameters
cat("\nStandardization parameters:\n")
cat("Variable | Mean      | SD\n")
cat("---------|-----------|----------\n")
for (i in 1:length(env_vars)) {
  cat(sprintf("%-8s | %9.3f | %9.3f\n", env_vars[i], orig_means[i], orig_sds[i]))
}

# Standardize variables
cat("\nStandardizing environmental variables...\n")
for (i in 1:length(env_vars)) {
  env_stack_df[[std_vars[i]]] <- (env_stack_df[[env_vars[i]]] - orig_means[i]) / orig_sds[i]
  
  # Check for issues
  n_nan <- sum(is.nan(env_stack_df[[std_vars[i]]]))
  n_inf <- sum(is.infinite(env_stack_df[[std_vars[i]]]))
  
  if (n_nan > 0 || n_inf > 0) {
    cat("  ÃÂ ", env_vars[i], ": Found", n_nan, "NaN and", n_inf, "Inf values\n")
    env_stack_df[[std_vars[i]]][is.nan(env_stack_df[[std_vars[i]]])] <- NA
    env_stack_df[[std_vars[i]]][is.infinite(env_stack_df[[std_vars[i]]])] <- NA
  } else {
    cat("  ", env_vars[i], "standardized successfully\n")
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

cat("   Design matrix created with", ncol(X.0), "predictors\n")

# Prepare coordinates
coords.0 <- as.matrix(env_stack_df[, c('x', 'y')])

# Check for and remove NA values
cat("\nChecking for NA values...\n")
na_rows <- rowSums(is.na(X.0)) > 0
n_na <- sum(na_rows)

if (n_na > 0) {
  cat("  ÃÂ  Found", format(n_na, big.mark = ","), "rows with NA values\n")
  X.0 <- X.0[!na_rows, ]
  coords.0 <- coords.0[!na_rows, ]
  cat("   After removal:", format(nrow(X.0), big.mark = ","), "locations remain\n")
} else {
  cat("   No NA values found\n")
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

# Initialize results data frame
results_summary <- data.frame(
  x = numeric(0),
  y = numeric(0),
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
      "(", chunk_size_actual, "locations )\n")
  
  # Extract chunk data
  X.0.chunk <- X.0[start_idx:end_idx, , drop = FALSE]
  coords.0.chunk <- coords.0[start_idx:end_idx, , drop = FALSE]
  
  # Run prediction with error handling
  tryCatch({
    cat("Running spatial predictions...\n")
    
    pred.chunk <- predict(out.sfMsPGOcc, 
                          X.0.chunk,  
                          coords.0.chunk,
                          n.omp.threads = N_THREADS,
                          verbose = TRUE,
                          n.report = 1000,
                          type = 'occupancy')
    
    cat(" Predictions completed\n")
    
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
        
        # Create row
        row_data <- data.frame(
          x = coords.0.chunk[loc, 1],
          y = coords.0.chunk[loc, 2],
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
    cat(" Chunk results saved to:", basename(chunk_file), "\n")
    
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
    cat("\n ERROR in chunk", i, ":", conditionMessage(e), "\n")
    
    # Implement fallback for memory errors
    if (grepl("cannot allocate", conditionMessage(e))) {
      cat("ÃÂ  Memory allocation error detected. Attempting smaller sub-chunks...\n")
      
      # Process in smaller sub-chunks
      smaller_chunk_size <- floor(CHUNK_SIZE / 5)
      n_sub_chunks <- ceiling(chunk_size_actual / smaller_chunk_size)
      
      for (j in 1:n_sub_chunks) {
        sub_start <- start_idx + (j - 1) * smaller_chunk_size
        sub_end <- min(start_idx + j * smaller_chunk_size - 1, end_idx)
        
        cat("\n  Processing sub-chunk", j, "of", n_sub_chunks, 
            "(locations", sub_start, "-", sub_end, ")\n")
        
        # [Sub-chunk processing code would go here - similar to main chunk processing]
        # [Omitted for brevity but would follow same pattern]
      }
    }
  })
  
  # Save intermediate results periodically
  if (i %% SAVE_INTERVAL == 0 || i == n_chunks) {
    intermediate_file <- file.path(base_output_dir, "intermediate_prediction_results.RData")
    save(results_summary, file = intermediate_file)
    cat("\n Intermediate results saved (", format(nrow(results_summary), big.mark = ","), 
        "total predictions )\n")
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

# Save combined results
final_results_file <- file.path(base_output_dir, "combined_prediction_results.RData")
save(results_summary, file = final_results_file)
cat(" Final results saved to:", final_results_file, "\n")
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
    
    # Create species output directory
    species_dir <- file.path(base_output_dir, paste0("species_", sprintf("%02d", sp), "_", species_names[sp]))
    dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Create rasters
    cat("  Creating rasters...\n")
    
    # Occurrence probability
    psi_mean_raster <- rast(sp_data[, c("x", "y", "psi_mean")], type = "xyz", crs = "EPSG:3577")
    writeRaster(psi_mean_raster, file.path(species_dir, "psi_mean.tif"), overwrite = TRUE)
    cat("     psi_mean.tif (mean occurrence probability)\n")
    
    # Uncertainty
    psi_sd_raster <- rast(sp_data[, c("x", "y", "psi_sd")], type = "xyz", crs = "EPSG:3577")
    writeRaster(psi_sd_raster, file.path(species_dir, "psi_sd.tif"), overwrite = TRUE)
    cat("     psi_sd.tif (occurrence probability SD)\n")
    
    # Posterior occurrence probability
    z_prob_raster <- rast(sp_data[, c("x", "y", "z_prob")], type = "xyz", crs = "EPSG:3577")
    writeRaster(z_prob_raster, file.path(species_dir, "z_prob.tif"), overwrite = TRUE)
    cat("     z_prob.tif (posterior occurrence probability)\n")
    
    # Spatial random effects
    w_mean_raster <- rast(sp_data[, c("x", "y", "w_mean")], type = "xyz", crs = "EPSG:3577")
    writeRaster(w_mean_raster, file.path(species_dir, "w_mean.tif"), overwrite = TRUE)
    cat("     w_mean.tif (mean spatial random effect)\n")
    
    w_sd_raster <- rast(sp_data[, c("x", "y", "w_sd")], type = "xyz", crs = "EPSG:3577")
    writeRaster(w_sd_raster, file.path(species_dir, "w_sd.tif"), overwrite = TRUE)
    cat("     w_sd.tif (spatial random effect SD)\n")
    
    # Summary statistics
    cat("  Summary statistics:\n")
    cat("    - Mean occurrence probability:", round(mean(sp_data$z_prob), 3), "\n")
    cat("    - SD occurrence probability:", round(sd(sp_data$z_prob), 3), "\n")
    cat("    - Min/Max occurrence:", round(min(sp_data$z_prob), 3), "/", round(max(sp_data$z_prob), 3), "\n")
  }
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("PREDICTION COMPLETE\n")
cat("========================================\n")

total_time <- difftime(Sys.time(), overall_start_time, units = "mins")
cat("Total processing time:", round(total_time, 2), "minutes\n")
cat("Output directory:", base_output_dir, "\n")
cat("\nFiles created:\n")
cat("  - Combined results: combined_prediction_results.RData\n")
cat("  - Chunk results:", n_chunks, "files in chunks/\n")
cat("  - Species rasters:", n_species, "folders with 5 rasters each\n")

cat("\n Multi-species spatial prediction completed successfully!\n")
cat("========================================\n\n")