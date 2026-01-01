##############################
#
#   Figure 3: Stopover Density on Drought Maps (2018)
#   Ian Becker
#   December 2025
#
##############################

# This script is used to make figure 3 which shows stopover
# by radar between drought conditions for fall and spring 2018

library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(ggspatial)
library(tigris)
library(ggnewscale)

options(tigris_use_cache = TRUE)

#################################
# DATA PREP
#################################

# Year for analysis

analysis_year <- 2018

# File paths 

fall_stopover_file <- "raster_data/stopover_fall.tif"
spring_stopover_file <- "raster_data/stopover_spring.tif"
fall_spei_file <- "SPEI_fall_1995_2020_250m.tif"  
spring_spei_file <- "raster_data/cropped_rasters/spring/SPEI_spring_1995_2020_spring_cropped.tif"
bcr_gdb <- "NABCI_ecoregion.gdb"

# Output settings (can change if there is a designated figure folder)

output_dir <- getwd()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Visualization parameters

platform_radius_km <- 80         # Radar detection radius (80 km)
platform_elevation_deg <- 0.4    # Visual elevation of platforms
platform_shadow_alpha <- 0.3     # Shadow transparency

# Plot settings

plot_width <- 8
plot_height <- 7
plot_dpi <- 600

#################################
# DATA LOAD
#################################

cat("Loading spatial data...\n")

# Load rasters

fall_stopover <- rast(fall_stopover_file)
spring_stopover <- rast(spring_stopover_file)
fall_spei <- rast(fall_spei_file)
spring_spei <- rast(spring_spei_file)

# Load boundaries

bcr <- st_read(bcr_gdb, layer = "BCR_Terrestrial_Master", quiet = TRUE)
states <- states(cb = TRUE, resolution = "20m")
tx_la <- states %>% filter(STUSPS %in% c("TX", "LA"))

# Transform to WGS84

flat_crs <- st_crs(4326)
bcr <- st_transform(bcr, crs(fall_stopover))
study_extent <- ext(fall_stopover)
bcr_crop <- st_crop(bcr, study_extent)
bcr_crop <- st_transform(bcr_crop, flat_crs)
tx_la <- st_transform(tx_la, flat_crs)
bcr_clipped <- st_intersection(bcr_crop, tx_la)

cat("Data loaded\n\n")

#################################
# RADAR STATION LOCATIONS
#################################

# Texas radars (n=12)

texas_radars <- data.frame(
  station = c("KGRK", "KHGX", "KCRP", "KDFX", "KFWS", "KDYX", 
              "KEPZ", "KLBB", "KMAF", "KSJT", "KAMA", "KBRO"),
  lat = c(30.722, 29.472, 27.784, 29.273, 32.573, 32.538,
          31.873, 33.654, 31.943, 31.371, 35.233, 25.916),
  lon = c(-97.383, -95.079, -97.511, -100.280, -97.303, -99.254,
          -106.698, -101.814, -102.189, -100.492, -101.709, -97.419)
)

# Louisiana radars (n=3)

louisiana_radars <- data.frame(
  station = c("KLCH", "KLIX", "KSHV"),
  lat = c(30.125, 30.337, 32.451),
  lon = c(-93.216, -89.826, -93.841)
)

# Combine (n=15 total)

radar_stations <- bind_rows(texas_radars, louisiana_radars)

# Check radars loaded correctly

cat("Radar stations:", nrow(radar_stations), "\n")
cat("  Texas:", nrow(texas_radars), "\n")
cat("  Louisiana:", nrow(louisiana_radars), "\n\n")

#################################
# FUNCTION: CREATE PLATFORM GEOMETRY
#################################

create_platform_geometry <- function(center_lon, center_lat, radius_km, elevation_deg) {
  
  # Convert km to degrees (approximate)
  
  radius_deg <- radius_km / 111
  n_points <- 100
  
  # Circle angles
  
  angles <- seq(0, 2*pi, length.out = n_points)
  
  # Platform shadow (on ground, slightly offset)
  
  shadow <- data.frame(
    x = center_lon + radius_deg * cos(angles) + 0.08,
    y = center_lat + radius_deg * sin(angles) - 0.08,
    component = "shadow",
    station = paste(center_lon, center_lat, sep="_")
  )
  
  # Vertical support line from ground to platform center
  
  support_line <- data.frame(
    x = c(center_lon, center_lon),
    y = c(center_lat, center_lat + elevation_deg),
    component = "support",
    station = paste(center_lon, center_lat, sep="_")
  )
  
  # Platform outline (elevated circle)
  
  outline <- data.frame(
    x = center_lon + radius_deg * cos(angles),
    y = center_lat + radius_deg * sin(angles) + elevation_deg,
    component = "outline",
    station = paste(center_lon, center_lat, sep="_")
  )
  
  return(list(shadow = shadow, support = support_line, outline = outline))
}

#################################
# FUNCTION: EXTRACT STOPOVER RASTER
#################################

extract_platform_stopover <- function(stopover_raster, center_lon, center_lat, 
                                      radius_km, elevation_deg) {
  
  # Create circular buffer
  
  center_sf <- st_sfc(st_point(c(center_lon, center_lat)), crs = 4326)
  buffer_sf <- st_buffer(center_sf, dist = radius_km * 1000)  # km to m
  
  # Crop and mask stopover to circle
  
  stopover_crop <- crop(stopover_raster, vect(buffer_sf))
  stopover_masked <- mask(stopover_crop, vect(buffer_sf))
  
  # Fill NA with 0 to show complete platform
  
  values(stopover_masked)[is.na(values(stopover_masked))] <- 0
  
  # Convert to dataframe
  
  stopover_df <- as.data.frame(stopover_masked, xy = TRUE)
  names(stopover_df)[3] <- "stopover"
  
  # Elevate to platform height
  
  stopover_df$y <- stopover_df$y + elevation_deg
  stopover_df$station <- paste(center_lon, center_lat, sep="_")
  
  return(stopover_df)
}

#################################
# FUNCTION: EXTRACT YEAR DATA
#################################

extract_year_data <- function(stopover_stack, spei_stack, year, season) {
  
  cat("  Extracting", season, year, "data...\n")
  
  # Find year in layer names
  
  stopover_years <- as.numeric(gsub(".*_(\\d{4})_.*", "\\1", names(stopover_stack)))
  spei_years <- as.numeric(gsub(".*_(\\d{4})_.*", "\\1", names(spei_stack)))
  
  stopover_idx <- which(stopover_years == year)
  spei_idx <- which(spei_years == year)
  
  if(length(stopover_idx) == 0 || length(spei_idx) == 0) {
    stop("Year ", year, " not found in ", season, " data")
  }
  
  # Extract and project to WGS84
  
  stopover_layer <- project(stopover_stack[[stopover_idx]], "EPSG:4326")
  spei_layer <- project(spei_stack[[spei_idx]], "EPSG:4326")
  
  return(list(stopover = stopover_layer, spei = spei_layer))
}

#################################
# FUNCTION: CREATE PANEL PLOT
#################################

create_panel <- function(stopover_raster, spei_raster, panel_label, season) {
  
  cat("  Building", season, "panel...\n")
  
  # Create platform geometry for all radars
  
  all_shadows <- data.frame()
  all_supports <- data.frame()
  all_outlines <- data.frame()
  
  for(i in 1:nrow(radar_stations)) {
    geom <- create_platform_geometry(
      center_lon = radar_stations$lon[i],
      center_lat = radar_stations$lat[i],
      radius_km = platform_radius_km,
      elevation_deg = platform_elevation_deg
    )
    all_shadows <- bind_rows(all_shadows, geom$shadow)
    all_supports <- bind_rows(all_supports, geom$support)
    all_outlines <- bind_rows(all_outlines, geom$outline)
  }
  
  # Extract stopover data on platforms
  
  all_stopover <- data.frame()
  for(i in 1:nrow(radar_stations)) {
    stopover_df <- extract_platform_stopover(
      stopover_raster,
      center_lon = radar_stations$lon[i],
      center_lat = radar_stations$lat[i],
      radius_km = platform_radius_km,
      elevation_deg = platform_elevation_deg
    )
    all_stopover <- bind_rows(all_stopover, stopover_df)
  }
  
  # Prepare SPEI background
  
  tx_la_buffered <- st_buffer(tx_la, dist = 0.5)
  spei_masked <- mask(spei_raster, vect(tx_la_buffered))
  spei_df <- as.data.frame(spei_masked, xy = TRUE, na.rm = TRUE)
  names(spei_df)[3] <- "spei"
  
  # Build plot
  
  p <- ggplot() +
    
    # Layer 1: SPEI drought background
    
    geom_raster(data = spei_df, aes(x = x, y = y, fill = spei)) +
    scale_fill_gradientn(
      colors = c("#8B4513", "#CD853F", "#DEB887", "#F5DEB3", 
                 "#F0F8FF", "#B0E0E6", "#4682B4", "#1E90FF"),
      name = "SPEI\n(Drought Index)",
      limits = c(-2.5, 2.5),
      breaks = c(-2, -1, 0, 1, 2),
      labels = c("-2 Extreme Drought", "-1 Dry", "0 Normal", "1 Wet", "2 Very Wet")
    ) +
    
    # Layer 2: State boundaries
    
    geom_sf(data = tx_la, fill = NA, color = "gray20", linewidth = 0.8) +
    
    # Layer 3: BCR boundaries
    
    geom_sf(data = bcr_clipped, fill = NA, color = "gray50", 
            linewidth = 0.25, linetype = "dotted", alpha = 0.6) +
    
    # Layer 4: Platform shadows
    
    geom_polygon(data = all_shadows,
                 aes(x = x, y = y, group = station),
                 fill = "black", color = NA, alpha = platform_shadow_alpha) +
    
    # Layer 5: Vertical support lines
    
    geom_line(data = all_supports,
              aes(x = x, y = y, group = station),
              color = "black", linewidth = 0.5, alpha = 0.7) +
    
    # Layer 6: Stopover raster on platforms (new fill scale)
    
    ggnewscale::new_scale_fill() +
    geom_raster(data = all_stopover, aes(x = x, y = y, fill = stopover)) +
    scale_fill_gradientn(
      colors = c("#fde724", "#7ad151", "#22a884", "#2a788e", "#414487", "#440154"),
      name = "Stopover\nDensity",
      na.value = "transparent"
    ) +
    
    # Layer 7: Platform outlines (solid black circles)
    
    geom_path(data = all_outlines,
              aes(x = x, y = y, group = station),
              color = "black", linewidth = 1.2) +
    
    # Annotations
    
    annotation_north_arrow(
      location = "tr", 
      which_north = "true",
      style = north_arrow_fancy_orienteering,
      height = unit(1.2, "cm"),
      width = unit(1.2, "cm")
    ) +
    annotation_scale(location = "bl", width_hint = 0.25, text_cex = 0.8) +
    
    # Labels and theme
    
    labs(
      title = panel_label,
      x = "Longitude",
      y = "Latitude"
    ) +
    coord_sf(xlim = c(-108.5, -88.5), ylim = c(24, 37), 
             crs = flat_crs, expand = FALSE) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      legend.position = "right",
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

#################################
# GENERATE PLOTS
#################################

# Panel (a): Fall 2018 (Wet conditions)

cat("Panel (a): Fall", analysis_year, "\n")
fall_data <- extract_year_data(fall_stopover, fall_spei, analysis_year, "fall")
p_fall <- create_panel(
  stopover_raster = fall_data$stopover,
  spei_raster = fall_data$spei,
  panel_label = paste0("(a) Fall ", analysis_year),
  season = "fall"
)

fall_file <- file.path(output_dir, paste0("panel_a_fall_", analysis_year, ".png"))
ggsave(fall_file, p_fall, width = plot_width, height = plot_height, dpi = plot_dpi)
cat("  ✓ Saved:", fall_file, "\n\n")

# Panel (b): Spring 2018 (Drought conditions)

cat("Panel (b): Spring", analysis_year, "\n")
spring_data <- extract_year_data(spring_stopover, spring_spei, analysis_year, "spring")
p_spring <- create_panel(
  stopover_raster = spring_data$stopover,
  spei_raster = spring_data$spei,
  panel_label = paste0("(b) Spring ", analysis_year),
  season = "spring"
)

spring_file <- file.path(output_dir, paste0("panel_b_spring_", analysis_year, ".png"))
ggsave(spring_file, p_spring, width = plot_width, height = plot_height, dpi = plot_dpi)
cat("Saved:", spring_file, "\n\n")

cat("FIGURE 3 COMPLETE!\n")
