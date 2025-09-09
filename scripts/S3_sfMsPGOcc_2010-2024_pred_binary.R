#!/usr/bin/env Rscript
# ============================================================================
# Script: Convert Continuous Habitat Predictions to Binary Maps using TSS
# Author: Adapted for Lamuel C.H. Chung's workflow
# Date: 2025-01-10
# Description: Converts continuous habitat suitability predictions to binary
#              habitat/non-habitat maps using TSS threshold.
#              Designed to work with climate projection outputs from sfMsPGOcc.
# ============================================================================

# Load required libraries
cat("\n========================================\n")
cat("LOADING REQUIRED LIBRARIES\n")
cat("========================================\n")
library(terra)
library(dplyr)
library(PresenceAbsence) # For TSS calculation
library(pROC)          # Alternative for ROC/TSS

# ============================================================================
# CONFIGURATION - Aligned with your existing structure
# ============================================================================
cat("\n========================================\n")
cat("CONFIGURATION\n")
cat("========================================\n")

# Location settings - matching your existing structure
LOCATION <- "mainland"     # Options: "mainland", "tasmania"
MODEL_TYPE <- "sfMsPGOcc"  # Your model type
PRED_FOLDER <- "pred_20250825" # The specific prediction run folder

# Base directories - matching your path structure
INPUT_BASE_DIR <- paste0("output/", LOCATION, "/", MODEL_TYPE, "/predictions")
OUTPUT_BASE_DIR <- paste0("output/", LOCATION, "/", MODEL_TYPE, "/binary")

# Species information - from your scripts
SPECIES_NAMES <- c("Dasyurus") # Spotted-tailed quoll
SPECIES_ID <- 2                 # Adjust based on which species from the multi-species model

# Climate scenarios - matching your script structure
MODELS <- c("gfdl-esm4", "mpi-esm1-2-hr", "ukesm1-0-ll", "mri-esm2-0", "ipsl-cm6a-lr")
SSPS <- c("ssp585", "ssp126", "ssp370")
TIME_PERIODS <- c("2071-2100", "2011-2040", "2041-2070")

# Create output directory
dir.create(OUTPUT_BASE_DIR, recursive = TRUE, showWarnings = FALSE)
cat("✓ Output directory created:", OUTPUT_BASE_DIR, "\n")

# ============================================================================
# LOAD OCCURRENCE DATA
# ============================================================================
cat("\n========================================\n")
cat("LOADING OCCURRENCE DATA\n")
cat("========================================\n")

# Load the original occurrence data used for model fitting
load("input/mainland/data_top12_2010-2024_5r.RData")

# Extract occurrence data for the target species
if (exists("data")) {
  # Get site coordinates
  site_coords <- data$coords
  
  # Get detection/non-detection data for the species
  y_data <- data$y[SPECIES_ID, , ] # [species, sites, replicates]
  
  # Create presence/absence at site level (any detection across replicates)
  site_presence <- apply(y_data, 1, function(x) any(x == 1, na.rm = TRUE))
  
  # Separate presence and absence coordinates
  presence_coords <- site_coords[site_presence, ]
  absence_coords <- site_coords[!site_presence, ]
  
  cat("✓ Occurrence data loaded\n")
  cat("  Presence sites:", nrow(presence_coords), "\n")
  cat("  Absence sites:", nrow(absence_coords), "\n")
} else {
  cat("⚠ Warning: Could not load occurrence data. Will use fixed threshold.\n")
  presence_coords <- NULL
  absence_coords <- NULL
}

# ============================================================================
# FUNCTIONS
# ============================================================================

#' Calculate optimal TSS threshold
#'
#' @param pred_raster SpatRaster of predicted probabilities
#' @param presence_coords Matrix of presence coordinates
#' @param absence_coords Matrix of absence coordinates
#' @param method Method for threshold selection ("MaxTSS", "MaxKappa", "MinDiff")
#' @return List with optimal threshold and performance metrics
calculate_optimal_threshold <- function(pred_raster, presence_coords, 
                                        absence_coords, method = "MaxTSS") {
  
  cat("Calculating optimal threshold using method:", method, "\n")
  
  # Extract values at presence and absence locations
  pres_vals <- terra::extract(pred_raster, presence_coords)[[1]]
  abs_vals <- terra::extract(pred_raster, absence_coords)[[1]]
  
  # Remove NA values
  pres_vals <- pres_vals[!is.na(pres_vals)]
  abs_vals <- abs_vals[!is.na(abs_vals)]
  
  cat("  Extracted values - Presences:", length(pres_vals), 
      "Absences:", length(abs_vals), "\n")
  
  # Create data for threshold optimization
  obs <- c(rep(1, length(pres_vals)), rep(0, length(abs_vals)))
  pred <- c(pres_vals, abs_vals)
  
  # Calculate performance across thresholds
  thresholds <- seq(min(pred, na.rm=T), max(pred, na.rm=T), length.out = 100)
  performance <- data.frame(
    threshold = thresholds,
    sensitivity = NA,
    specificity = NA,
    tss = NA,
    kappa = NA
  )
  
  for (i in seq_along(thresholds)) {
    thresh <- thresholds[i]
    pred_binary <- as.numeric(pred >= thresh)
    
    # Confusion matrix
    tp <- sum(pred_binary == 1 & obs == 1, na.rm = TRUE)
    tn <- sum(pred_binary == 0 & obs == 0, na.rm = TRUE)
    fp <- sum(pred_binary == 1 & obs == 0, na.rm = TRUE)
    fn <- sum(pred_binary == 0 & obs == 1, na.rm = TRUE)
    
    # Calculate metrics
    if((tp + fn) > 0) performance$sensitivity[i] <- tp / (tp + fn)
    if((tn + fp) > 0) performance$specificity[i] <- tn / (tn + fp)
    if(!is.na(performance$sensitivity[i]) && !is.na(performance$specificity[i])) {
      performance$tss[i] <- performance$sensitivity[i] + performance$specificity[i] - 1
    }
    
    # Kappa
    po <- (tp + tn) / length(obs)
    pe <- ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / length(obs)^2
    if((1-pe) != 0) performance$kappa[i] <- (po - pe) / (1 - pe)
  }
  
  # Select optimal threshold based on method
  if (method == "MaxTSS") {
    optimal_idx <- which.max(performance$tss)
  } else if (method == "MaxKappa") {
    optimal_idx <- which.max(performance$kappa)
  } else if (method == "MinDiff") {
    diff <- abs(performance$sensitivity - performance$specificity)
    optimal_idx <- which.min(diff)
  }
  
  optimal <- performance[optimal_idx, ]
  
  cat("  Optimal threshold:", round(optimal$threshold, 4), "\n")
  cat("  TSS:", round(optimal$tss, 3), "\n")
  cat("  Sensitivity:", round(optimal$sensitivity, 3), "\n")
  cat("  Specificity:", round(optimal$specificity, 3), "\n")
  
  return(list(
    threshold = optimal$threshold,
    performance = optimal,
    all_performance = performance
  ))
}


#' Process predictions for a specific scenario by building the exact file path.
#'
#' @param base_scenario_dir Directory for the climate scenario (e.g., .../baseline_1981-2010)
#' @param output_dir Output directory for binary rasters
#' @param threshold Threshold for binary conversion
#' @param species_id Numeric ID of the species
#' @param species_name Character name of the species (e.g., "Dasyurus")
#' @param scenario_name Character name of the scenario (e.g., "baseline_1981-2010")
#' @param pred_folder The specific prediction run folder name (e.g., "pred_20250805")
process_scenario <- function(base_scenario_dir, output_dir, threshold, species_id, 
                             species_name, scenario_name, pred_folder) {
  
  # 1. Construct the expected folder and file path based on your structure
  species_folder <- paste0("species_", sprintf("%02d", species_id), "_", species_name, "_", scenario_name)
  file_name <- paste0("z_prob_", scenario_name, ".tif")
  
  # This is the full, direct path to the target raster
  target_file <- file.path(base_scenario_dir, pred_folder, species_folder, file_name)
  
  cat("  -> Targeting file:", target_file, "\n")
  
  # 2. Check if the specific file exists
  if (!file.exists(target_file)) {
    cat("    ⚠ File not found. Skipping.\n")
    return(NULL)
  }
  
  cat("    ✓ File found. Processing...\n")
  
  # 3. Process the file
  # Read prediction raster
  pred_raster <- terra::rast(target_file)
  
  # Convert to binary
  binary_raster <- pred_raster >= threshold
  
  # Create output filename
  base_name <- tools::file_path_sans_ext(basename(target_file))
  output_file <- file.path(output_dir, paste0(base_name, "_binary_tss", 
                                              round(threshold * 100), ".tif"))
  
  # Save binary raster
  terra::writeRaster(binary_raster, output_file, overwrite = TRUE)
  cat("      ✓ Saved:", basename(output_file), "\n")
  
  return(output_dir)
}

# ============================================================================
# MAIN PROCESSING
# ============================================================================
cat("\n========================================\n")
cat("STARTING BINARY CONVERSION\n")
cat("========================================\n")

# Step 1: Calculate optimal threshold using baseline predictions
cat("\nStep 1: Calculating optimal TSS threshold\n")
cat("----------------------------------------\n")

optimal_threshold <- 0.5 # Default value

if (!is.null(presence_coords) && !is.null(absence_coords)) {
  # Dynamically construct the path to the baseline prediction file
  baseline_scenario_name <- "baseline_1981-2010"
  species_folder <- paste0("species_", sprintf("%02d", SPECIES_ID), "_", SPECIES_NAMES[1], "_", baseline_scenario_name)
  file_name <- paste0("z_prob_", baseline_scenario_name, ".tif")
  
  baseline_pred_file <- file.path(INPUT_BASE_DIR, "climate_projections", baseline_scenario_name, PRED_FOLDER, species_folder, file_name)
  
  cat("Looking for baseline prediction file for thresholding:\n ->", baseline_pred_file, "\n")
  
  if (file.exists(baseline_pred_file)) {
    sample_raster <- terra::rast(baseline_pred_file)
    
    # Calculate optimal threshold
    threshold_result <- calculate_optimal_threshold(
      sample_raster, 
      presence_coords, 
      absence_coords,
      method = "MaxTSS"
    )
    
    optimal_threshold <- threshold_result$threshold
    
    # Save threshold information
    threshold_info <- data.frame(
      location = LOCATION,
      model = MODEL_TYPE,
      species = SPECIES_NAMES[1],
      threshold = optimal_threshold,
      tss = threshold_result$performance$tss,
      sensitivity = threshold_result$performance$sensitivity,
      specificity = threshold_result$performance$specificity,
      date_processed = Sys.Date()
    )
    
    write.csv(threshold_info, 
              file.path(OUTPUT_BASE_DIR, "threshold_info.csv"), 
              row.names = FALSE)
    
    # Create performance plot
    pdf(file.path(OUTPUT_BASE_DIR, "threshold_optimization.pdf"), 
        width = 10, height = 8)
    par(mfrow = c(2, 2))
    
    # TSS plot
    plot(threshold_result$all_performance$threshold, 
         threshold_result$all_performance$tss,
         type = "l", lwd = 2, col = "blue",
         xlab = "Threshold", ylab = "TSS",
         main = "True Skill Statistic")
    abline(v = optimal_threshold, col = "red", lty = 2)
    
    # Sensitivity and Specificity
    plot(threshold_result$all_performance$threshold, 
         threshold_result$all_performance$sensitivity,
         type = "l", lwd = 2, col = "darkgreen",
         xlab = "Threshold", ylab = "Rate",
         main = "Sensitivity vs Specificity", ylim = c(0, 1))
    lines(threshold_result$all_performance$threshold, 
          threshold_result$all_performance$specificity,
          lwd = 2, col = "darkorange")
    abline(v = optimal_threshold, col = "red", lty = 2)
    legend("right", c("Sensitivity", "Specificity"), 
           col = c("darkgreen", "darkorange"), lwd = 2)
    
    # Kappa
    plot(threshold_result$all_performance$threshold, 
         threshold_result$all_performance$kappa,
         type = "l", lwd = 2, col = "purple",
         xlab = "Threshold", ylab = "Kappa",
         main = "Cohen's Kappa")
    abline(v = optimal_threshold, col = "red", lty = 2)
    
    dev.off()
    cat("✓ Threshold optimization complete. Plot saved.\n")
    
  } else {
    cat("⚠ No baseline prediction file found at the specified path. Using default threshold.\n")
  }
} else {
  cat("⚠ No occurrence data available. Using default threshold.\n")
}

cat("\n✓ Using threshold:", round(optimal_threshold, 4), "\n")

# Step 2: Process all scenarios
cat("\nStep 2: Processing all climate scenarios\n")
cat("----------------------------------------\n")

# Process baseline
cat("\nProcessing baseline (1981-2010)...\n")
baseline_scenario_name <- "baseline_1981-2010"
base_scenario_dir <- file.path(INPUT_BASE_DIR, "climate_projections", baseline_scenario_name)
if (dir.exists(base_scenario_dir)) {
  output_dir <- file.path(OUTPUT_BASE_DIR, baseline_scenario_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  process_scenario(base_scenario_dir, output_dir, optimal_threshold, 
                   SPECIES_ID, SPECIES_NAMES[1], baseline_scenario_name, PRED_FOLDER)
}

# Process future projections
for (period in TIME_PERIODS) {
  for (ssp in SSPS) {
    for (model in MODELS) {
      scenario_name <- paste(period, model, ssp, sep = "_")
      # --- ADD THIS LINE FOR DEBUGGING ---
      cat("--> Inspecting scenario_name: '", scenario_name, "'\n", sep = "")
      cat("\nProcessing:", scenario_name, "\n")
      
      
      
      # Base input directory for this scenario
      base_scenario_dir <- file.path(INPUT_BASE_DIR, "climate_projections", scenario_name)
      
      if (!dir.exists(base_scenario_dir)) {
        cat("  ⚠ Base directory not found:", base_scenario_dir, "\n")
        next
      }
      
      # Output directory for this scenario
      output_dir <- file.path(OUTPUT_BASE_DIR, scenario_name)
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Process the scenario using the updated function
      process_scenario(base_scenario_dir, output_dir, optimal_threshold, 
                       SPECIES_ID, SPECIES_NAMES[1], scenario_name, PRED_FOLDER)
    }
  }
}

# Step 3: Create summary
cat("\n========================================\n")
cat("PROCESSING COMPLETE\n")
cat("========================================\n")
cat("Binary maps saved to:", OUTPUT_BASE_DIR, "\n")
cat("Threshold used:", round(optimal_threshold, 4), "\n")
cat("Check threshold_info.csv and threshold_optimization.pdf for details\n")
cat("========================================\n\n")