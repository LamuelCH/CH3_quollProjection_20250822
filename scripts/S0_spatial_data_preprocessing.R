# library(raster)
# library(magrittr)
# library(terra)
# library(slga)
# library(ggplot2)
# library(sf)
# library(geodata)
# library(tidyr)

###########################################################################
# Create species range mask ---------------------------------------------------
library(terra)
library(dplyr)

options(mc.core = parallel::detectCores())

# Define study area
# SE mainalnd & Tasmania AOI
aoi_SE_3577 = ext(684147.3, 2225320.5, -4370000, -2556437.4)
aoi_tas_3577 = ext(1000000, 1430000, -4850000, -4370000)


# Create WGS84 AOI --------------------------------------------------------
  # SE mainland -------------------------------------------------------------
  # Transform to WGS84 to get the equivalent extent
  {aoi_SE_wgs84_exact = as.polygons(aoi_SE_3577, crs = "EPSG:3577") %>% 
    project(.,"EPSG:4326") %>% 
    ext(.)
  
  buffer_deg <- 1
  
  aoi_SE_wgs84 <- ext(aoi_SE_wgs84_exact[1] - buffer_deg,  # xmin
                   aoi_SE_wgs84_exact[2] + buffer_deg,  # xmax  
                   aoi_SE_wgs84_exact[3] - buffer_deg,  # ymin
                   aoi_SE_wgs84_exact[4] + buffer_deg)  # ymax
  remove(aoi_SE_wgs84_exact)}


  # Tasmania ----------------------------------------------------------------
  {aoi_tas_wgs84_exact = as.polygons(aoi_tas_3577, crs = "EPSG:3577") %>% 
    project(.,"EPSG:4326") %>% 
    ext(.)
  
  buffer_deg <- 1
  
  aoi_tas_wgs84 <- ext(aoi_tas_wgs84_exact[1] - buffer_deg,  # xmin
                      aoi_tas_wgs84_exact[2] + buffer_deg,  # xmax  
                      aoi_tas_wgs84_exact[3] - buffer_deg,  # ymin
                      aoi_tas_wgs84_exact[4] + buffer_deg)  # ymax
  remove(aoi_tas_wgs84_exact)}
  
# Create a raster mask with 1 km resolution
# Generate an empty raster with the defined extent and resolution
raster_SE_3577 = rast(extent = aoi_SE_3577, resolution = 1000, crs = "EPSG:3577")
raster_tas_3577 = rast(extent = aoi_tas_3577, resolution = 1000, crs = "EPSG:3577")

# Foliage Projective Cover ------------------------------------------------
# Foliage projective cover
# SE mainaland
fpc_SE = rast("input/raw/Foliage Projective Cover/Woody vegetation cover - Landsat, JRSRP, Australian coverage, 2000-2010/lztmre_aus_y20002011_dm7a2_d20050630_r500m.tif")
fpc_SE <- terra::project(fpc_SE, raster_SE_3577, method = "bilinear") 
fpc_SE = clamp((fpc_SE - 100) / (200 - 100), 0, 1)
names(fpc_SE) = "fpc"
plot(fpc_SE)
writeRaster(fpc_SE, filename = "input/mainland/mainland_env_fpc_EPSG3577.tif", overwrite = TRUE)

# Tasmania
fpc_tas = rast("input/raw/Foliage Projective Cover/Woody vegetation cover - Landsat, JRSRP, Australian coverage, 2000-2010/lztmre_aus_y20002011_dm7a2_d20050630_r500m.tif")
fpc_tas <- terra::project(fpc_tas, raster_tas_3577, method = "bilinear")
fpc_tas = clamp((fpc_tas - 100) / (200 - 100), 0, 1)
names(fpc_tas) = "fpc"
plot(fpc_tas)
writeRaster(fpc_tas, filename = "input/tasmania/tasmania_env_fpc_EPSG3577.tif", overwrite = TRUE)

# Create binary mask from FPC
binary_SE <- terra::ifel(is.na(fpc_SE), NA, 1)
binary_tas = terra::ifel(is.na(fpc_tas), NA, 1)

###########################################################################
# CHELSA bioclimatic data -------------------------------------------------------
# Future -------------------------------------------------
# Define models and SSPs
models <- c("gfdl-esm4", "mpi-esm1-2-hr", "ukesm1-0-ll", "mri-esm2-0", "ipsl-cm6a-lr")
ssps <- c("ssp585", "ssp126", "ssp370")

# Define time periods
time_periods <- c("2071-2100", "2011-2040", "2041-2070")

# Define base input directory (without specific time period)
base_input_dir <- "input/raw/CHELSA/climatologies/"
base_output_dir <- "input/mainland/"

# Optional: Set terra options for better memory management
terraOptions(memfrac = 0.8)  # Use up to 80% of available RAM
terraOptions(tempdir = tempdir())  # Ensure temp files go to appropriate directory

# Loop through each time period, model, and SSP
for (time_period in time_periods) {
  cat("\nProcessing time period:", time_period, "\n")
  
  for (ssp in ssps) {
    cat("  Processing SSP:", ssp, "\n")
    
    # Create output directory for this time period and SSP
    output_dir <- file.path(base_output_dir, time_period, ssp)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    for (model in models) {
      cat("    Processing model:", model, "\n")
      
      # Construct input directory path
      input_dir <- file.path(base_input_dir, time_period, model, ssp, "bio")
      
      # List all raster files in the directory
      bio_files <- list.files(input_dir, pattern = "^CHELSA_bio", full.names = TRUE)
      
      # Proceed only if there are raster files
      if (length(bio_files) > 0) {
        
        # Define output file path
        output_file <- file.path(output_dir, paste0("mainland_CHELSA_bio_", time_period, "_", model, "_", ssp, "_V.2.1_EPSG3577.tif"))
        
        force_reprocess <- FALSE  # Set to TRUE to reprocess existing files
        
        # Check if output file already exists
        if (file.exists(output_file) && !force_reprocess) {
          cat("    SKIPPING", model, "- file already exists:", basename(output_file), "\n")
          next  # Skip to next iteration
        }
        
        tryCatch({
          # Extract layer names from file names
          layer_names <- sub(".*CHELSA_([^_]+)_.*", "\\1", basename(bio_files))
          
          # Load rasters and stack them
          bio <- rast(bio_files)
          
          # Assign names to layers
          names(bio) <- layer_names
          
          # OPTIMIZATION: Crop to AOI first (in WGS84) before reprojection
          bio_cropped <- crop(bio, aoi_SE_wgs84)
          
          # Clear original bio from memory
          rm(bio)
          
          # Reproject the cropped raster to EPSG:3577
          bio_reprojected <- terra::project(bio_cropped, raster_SE_3577, method = "bilinear")
          
          # Clear cropped bio from memory
          rm(bio_cropped)
          
          # Crop to exact target extent to ensure all outputs have identical extents
          bio_final <- crop(bio_reprojected, aoi_SE_3577)
          
          # Clear reprojected bio from memory
          rm(bio_reprojected)
          
          # Apply mask after final cropping
          bio_final <- bio_final * binary_SE
          
          # Save the processed raster stack (output_file already defined above)
          writeRaster(bio_final, output_file, overwrite = TRUE)
          
          # Clear final bio from memory
          rm(bio_final)
          
          print(paste("Processed and saved:", output_file))
          
        }, error = function(e) {
          cat("    ERROR processing", model, ":", e$message, "\n")
        })
        
        # MEMORY CLEANUP after each model
        # Force garbage collection
        gc(verbose = FALSE)
        
        # Clear any temporary files created by terra
        tmpFiles(remove = TRUE)
        
        # Optional: Print memory usage for monitoring
        # cat("    Memory used:", format(object.size(ls(envir = .GlobalEnv)), units = "MB"), "\n")
        
      } else {
        warning(paste("No raster files found for", time_period, model, ssp))
      }
    }
    
    # Additional cleanup after each SSP
    gc(verbose = FALSE)
    tmpFiles(remove = TRUE)
  }
  
  # Major cleanup after each time period
  gc(verbose = FALSE)
  tmpFiles(remove = TRUE)
  
  cat("Completed time period:", time_period, "\n")
}

# Final cleanup
gc(verbose = FALSE)
tmpFiles(remove = TRUE)
cat("\nAll processing completed successfully!\n")


  # Tasmania ----------------------------------------------------------------
  # Define models and SSPs
  # Define models and SSPs
  models <- c("gfdl-esm4", "mpi-esm1-2-hr", "ukesm1-0-ll", "mri-esm2-0", "ipsl-cm6a-lr")
  ssps <- c("ssp585", "ssp126", "ssp370")
  
  # Define time periods
  time_periods <- c("2071-2100", "2011-2040", "2041-2070")
  
  # Define base input directory (without specific time period)
  base_input_dir <- "input/raw/CHELSA/climatologies/"
  base_output_dir <- "input/tasmania/"
  
  # Optional: Set terra options for better memory management
  terraOptions(memfrac = 0.8)  # Use up to 80% of available RAM
  terraOptions(tempdir = tempdir())  # Ensure temp files go to appropriate directory
  
  # Loop through each time period, model, and SSP
  for (time_period in time_periods) {
    cat("\nProcessing time period:", time_period, "\n")
    
    for (ssp in ssps) {
      cat("  Processing SSP:", ssp, "\n")
      
      # Create output directory for this time period and SSP
      output_dir <- file.path(base_output_dir, time_period, ssp)
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
      
      for (model in models) {
        cat("    Processing model:", model, "\n")
        
        # Construct input directory path
        input_dir <- file.path(base_input_dir, time_period, model, ssp, "bio")
        
        # List all raster files in the directory
        bio_files <- list.files(input_dir, pattern = "^CHELSA_bio", full.names = TRUE)
        
        # Proceed only if there are raster files
        if (length(bio_files) > 0) {
          
          # Define output file path
          output_file <- file.path(output_dir, paste0("tasmania_CHELSA_bio_", time_period, "_", model, "_", ssp, "_V.2.1_EPSG3577.tif"))
          
          force_reprocess <- FALSE  # Set to TRUE to reprocess existing files
          
          # Check if output file already exists
          if (file.exists(output_file) && !force_reprocess) {
            cat("    SKIPPING", model, "- file already exists:", basename(output_file), "\n")
            next  # Skip to next iteration
          }
          
          tryCatch({
            # Extract layer names from file names
            layer_names <- sub(".*CHELSA_([^_]+)_.*", "\\1", basename(bio_files))
            
            # Load rasters and stack them
            bio <- rast(bio_files)
            
            # Assign names to layers
            names(bio) <- layer_names
            
            # OPTIMIZATION: Crop to AOI first (in WGS84) before reprojection
            bio_cropped <- crop(bio, aoi_tas_wgs84)
            
            # Clear original bio from memory
            rm(bio)
            
            # Reproject the cropped raster to EPSG:3577
            bio_reprojected <- terra::project(bio_cropped, raster_tas_3577, method = "bilinear")
            
            # Clear cropped bio from memory
            rm(bio_cropped)
            
            # Crop to exact target extent to ensure all outputs have identical extents
            bio_final <- crop(bio_reprojected, aoi_tas_3577)
            
            # Clear reprojected bio from memory
            rm(bio_reprojected)
            
            # Apply mask after final cropping
            bio_final <- bio_final * binary_tas
            
            # Save the processed raster stack (output_file already defined above)
            writeRaster(bio_final, output_file, overwrite = TRUE)
            
            # Clear final bio from memory
            rm(bio_final)
            
            print(paste("Processed and saved:", output_file))
            
          }, error = function(e) {
            cat("    ERROR processing", model, ":", e$message, "\n")
          })
          
          # MEMORY CLEANUP after each model
          # Force garbage collection
          gc(verbose = FALSE)
          
          # Clear any temporary files created by terra
          tmpFiles(remove = TRUE)
          
          # Optional: Print memory usage for monitoring
          # cat("    Memory used:", format(object.size(ls(envir = .GlobalEnv)), units = "MB"), "\n")
          
        } else {
          warning(paste("No raster files found for", time_period, model, ssp))
        }
      }
      
      # Additional cleanup after each SSP
      gc(verbose = FALSE)
      tmpFiles(remove = TRUE)
    }
    
    # Major cleanup after each time period
    gc(verbose = FALSE)
    tmpFiles(remove = TRUE)
    
    cat("Completed time period:", time_period, "\n")
  }
  
  # Final cleanup
  gc(verbose = FALSE)
  tmpFiles(remove = TRUE)
  cat("\nAll processing completed successfully!\n")


# Current -----------------------------------------------------------------
  # SE Mainland -------------------------------------------------------------
    # 1981-2010 CHELSA data
    ## Define the directory containing the CHELSA rasters
    {input_dir <- "input/raw/CHELSA/climatologies/1981-2010/bio"
    output_dir <- "input/mainland/1981-2010/"
    
    # List all raster files in the directory
    bio <- list.files(input_dir, pattern = "^CHELSA_bio", full.names = TRUE)
    
    # Extract layer names from file names
    layer_names <- sub(".*CHELSA_([^_]+)_.*", "\\1", basename(bio))
    
    bio = rast(bio)
    bio_cropped <- crop(bio, aoi_SE_wgs84) # pre-crop
    
    # Reproject the raster using bilinear method
    bio <- terra::project(bio_cropped, raster_SE_3577, method = "bilinear") * binary_SE
    
    # 4. Mask the raster.
    names(bio) <- layer_names
    
    # Save the stacked raster to a file
    output_file = file.path(output_dir, "mainland_CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
    writeRaster(bio,output_file, overwrite = TRUE)
    remove(bio_cropped)}
    
# Tasmania ----------------------------------------------------------------
    # 1981-2010 CHELSA data
    ## Define the directory containing the CHELSA rasters
    {input_dir <- "input/raw/CHELSA/climatologies/1981-2010/bio"
    output_dir <- "input/tasmania/1981-2010/"
    
    # List all raster files in the directory
    bio <- list.files(input_dir, pattern = "^CHELSA_bio", full.names = TRUE)
    
    # Extract layer names from file names
    layer_names <- sub(".*CHELSA_([^_]+)_.*", "\\1", basename(bio))
    
    bio = rast(bio)
    bio_cropped <- crop(bio, aoi_tas_wgs84) # pre-crop
    
    # Reproject the raster using bilinear method
    bio <- terra::project(bio_cropped, raster_tas_3577, method = "bilinear") * binary_tas
    
    # 4. Mask the raster.
    names(bio) <- layer_names
    
    # Save the stacked raster to a file
    output_file = file.path(output_dir, "tasmania_CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
    writeRaster(bio,output_file, overwrite = TRUE)
    remove(bio_cropped)}
  


###########################################################################
# Terrain layers -----------------------------------------------------------
  # SE mainland
  dem_SE = rast("input/raw/3secSRTM_DEM/DEM_ESRI_GRID_16bit_Integer/dem3s_int/hdr.adf") %>% crop(., aoi_SE_wgs84)
  names(dem_SE) = "dem"
  tpi_SE = terrain(dem_SE, "TPI")
  names(tpi_SE) = "tpi"
  terrain_SE = c(dem_SE, tpi_SE)
  terrain_SE = terra::project(terrain_SE, raster_SE_3577, method = "bilinear") * binary_SE
  
  writeRaster(terrain_SE, filename = "input/mainland/mainland_env_terrain_EPSG3577.tif", overwrite = TRUE)
  
  # Tasmania
  dem_tas = rast("input/raw/3secSRTM_DEM/DEM_ESRI_GRID_16bit_Integer/dem3s_int/hdr.adf") %>% crop(., aoi_tas_wgs84)
  names(dem_tas) = "dem"
  tpi_tas = terrain(dem_tas, "TPI")
  names(tpi_tas) = "tpi"
  terrain_tas = c(dem_tas, tpi_tas)
  terrain_tas = terra::project(terrain_tas, raster_tas_3577, method = "bilinear") * binary_tas
  
  writeRaster(terrain_tas, filename = "input/tasmania/tasmania_env_terrain_EPSG3577.tif", overwrite = TRUE)
  


# Test if all raster have same extent -------------------------------------
  bio_SE = rast("input/mainland/1981-2010/mainland_CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
  fpc_SE = rast("input/mainland/mainland_env_fpc_EPSG3577.tif")
  stack_SE = c(terrain_SE,fpc_SE,bio_SE )

  bio_tas = rast("input/tasmania/1981-2010/tasmania_CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
  fpc_tas = rast("input/tasmania/tasmania_env_fpc_EPSG3577.tif")
  stack_tas = c(terrain_tas,fpc_tas,bio_tas )
  

# Pearson Correlation Analysis --------------------------------------------
###########################################################################
library(ggplot2)
library(corrplot)

# Bioclimatic variables
bio = rast("input/raster/1981-2010/CHELSA_bio_1981-2010_V.2.1_EPSG3577.tif")
corr_bio = layerCor(bio_masked, "pearson", na.rm = TRUE)

bio = bio[c("bio5", "bio6")]

# terrain variables
terrain = rast("input/env_terrain_EPSG3577.tif")
corr_terrain = layerCor(terrain, "pearson", na.rm = TRUE)

# plot resutls 
corrplot(corr_terrain$correlation, 
         method = "circle", 
         type = "lower", 
         insig = "blank", 
         addCoef.col = "black", 
         number.cex = 0.6)

corrplot(corr_bio$correlation, 
         method = "circle", 
         type = "lower", 
         insig = "blank", 
         addCoef.col = "black", 
         number.cex = 0.6)


env_stack = c(subset(bio_masked, c("bio5", "bio6", "bio12", "bio15")), 
              subset(terrain, c("dem", "eastness", "northness", "tpi", "tri", "rock")))
correlation_matrix <- layerCor(env_stack, "pearson", na.rm = TRUE)
corrplot(correlation_matrix$correlation, method = "circle", type = "lower", insig = "blank", addCoef.col = "black", number.cex = 0.6)
