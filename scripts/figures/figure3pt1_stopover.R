##############################
#
#   Figure 3: Elevated Radar Platforms with Stopover
#   Ian Becker
#   12/18/2025
#
#   Creates elevated circular platforms at each radar location
#   with stopover patterns shown as texture on platform surfaces
#
##############################

library(terra)
library(sf)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggspatial)
library(tigris)
library(viridis)
library(ggnewscale)  # For multiple fill scales

options(tigris_use_cache = TRUE)
setwd("/Users/ianbecker/Library/CloudStorage/OneDrive-TheUniversityofTexas-RioGrandeValley/DroughtRadar")

#################################
# CONFIGURATION
#################################

# SELECT EXEMPLAR YEARS
fall_drought_year <- 2018
fall_wet_year <- 2018
spring_drought_year <- 2018
spring_wet_year <- 2010

# File paths
fall_stopover_file <- "raster_data/stopover_fall.tif"
spring_stopover_file <- "raster_data/stopover_spring.tif"
fall_spei_file <- "SPEI_fall_1995_2020_250m.tif"  
spring_spei_file <- "raster_data/cropped_rasters/spring_raster/SPEI_spring_1995_2020.tif"
bcr_gdb <- "NABCI_ecoregion.gdb"

# Output settings
output_dir <- "figure3_platforms"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Platform visualization settings
platform_radius_km <- 80              # Radar detection radius
platform_elevation_deg <- 0.4         # Visual elevation of platforms (REDUCED - closer to ground)
platform_side_color <- "#2c3e50"      # Dark color for platform sides
platform_shadow_alpha <- 0.3          # Shadow transparency
stopover_sample_distance <- 10000     # Sample every 10km for stopover spikes on platform

# Plot dimensions
plot_width <- 8
plot_height <- 7
plot_dpi <- 300

cat(paste(rep("=", 70), collapse=""), "\n")
cat("FIGURE 3: ELEVATED RADAR PLATFORMS WITH STOPOVER\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("Fall drought year:", fall_drought_year, "\n")
cat("Fall wet year:", fall_wet_year, "\n")
cat("Spring drought year:", spring_drought_year, "\n")
cat("Spring wet year:", spring_wet_year, "\n\n")

#################################
# LOAD BASE DATA
#################################

cat("Loading base spatial data...\n")

fall_stopover <- rast(fall_stopover_file)
spring_stopover <- rast(spring_stopover_file)
fall_spei <- rast(fall_spei_file)
spring_spei <- rast(spring_spei_file)

years <- 1995:(1995 + nlyr(fall_stopover) - 1)

bcr <- st_read(bcr_gdb, layer = "BCR_Terrestrial_Master", quiet = TRUE)
states <- states(cb = TRUE, resolution = "20m")
tx_la <- states %>% filter(STUSPS %in% c("TX", "LA"))

flat_crs <- st_crs(4326)

bcr <- st_transform(bcr, crs(fall_stopover))
study_extent <- ext(fall_stopover)
bcr_crop <- st_crop(bcr, study_extent)
bcr_crop <- st_transform(bcr_crop, flat_crs)
tx_la <- st_transform(tx_la, flat_crs)
bcr_clipped <- st_intersection(bcr_crop, tx_la)

cat("Base data loaded successfully!\n\n")

#################################
# RADAR STATION LOCATIONS
#################################

cat("Setting up radar station locations...\n")

texas_radars <- data.frame(
  site = c("KGRK", "KHGX", "KCRP", "KDFX", "KFWS", 
           "KDYX", "KEPZ", "KLBB", "KMAF", "KSJT", "KAMA", "KBRO"),
  name = c("Ft. Hood", "Houston", "Corpus Christi", 
           "Del Rio", "Dallas/Ft. Worth", "Dyess AFB", "El Paso", 
           "Lubbock", "Midland", "San Angelo", "Amarillo", "Brownsville"),
  latitude = c(30.722, 29.472, 27.784, 29.273, 32.573, 
               32.538, 31.873, 33.654, 31.943, 31.371, 35.233, 25.916),
  longitude = c(-97.383, -95.079, -97.511, -100.280, -97.303, 
                -99.254, -106.698, -101.814, -102.189, -100.492, -101.709, -97.419)
)

louisiana_radars <- data.frame(
  site = c("KLCH", "KLIX", "KSHV"),
  name = c("Lake Charles", "New Orleans", "Shreveport"),
  latitude = c(30.125, 30.337, 32.451),
  longitude = c(-93.216, -89.826, -93.841)
)

radar_stations <- bind_rows(texas_radars, louisiana_radars) %>%
  rename(station = site, lat = latitude, lon = longitude)

cat("  Radar stations:", nrow(radar_stations), "\n\n")

#################################
# FUNCTION: CREATE PLATFORM GEOMETRY
#################################

create_platform_geometry <- function(center_lon, center_lat, radius_km, elevation_deg, n_points = 100) {
  
  # Convert km to degrees (rough approximation)
  radius_deg <- radius_km / 111
  
  # Create angles for circle
  angles <- seq(0, 2*pi, length.out = n_points)
  
  # Platform shadow (on ground, slightly offset)
  shadow_offset_x <- 0.08
  shadow_offset_y <- -0.08
  shadow <- data.frame(
    x = center_lon + radius_deg * cos(angles) + shadow_offset_x,
    y = center_lat + radius_deg * sin(angles) + shadow_offset_y,
    component = "shadow",
    station = center_lon  # For grouping
  )
  
  # Platform top (elevated circle)
  top <- data.frame(
    x = center_lon + radius_deg * cos(angles),
    y = center_lat + radius_deg * sin(angles) + elevation_deg,
    component = "top",
    station = center_lon
  )
  
  # Platform sides - visible front edge (southern arc)
  front_angles <- angles[angles >= pi * 0.7 & angles <= pi * 1.3]
  n_front <- length(front_angles)
  
  side_front <- data.frame()
  for(i in 1:(n_front-1)) {
    # Create quad strip for this segment
    quad <- data.frame(
      x = c(center_lon + radius_deg * cos(front_angles[i]),
            center_lon + radius_deg * cos(front_angles[i+1]),
            center_lon + radius_deg * cos(front_angles[i+1]),
            center_lon + radius_deg * cos(front_angles[i])),
      y = c(center_lat + radius_deg * sin(front_angles[i]),
            center_lat + radius_deg * sin(front_angles[i+1]),
            center_lat + radius_deg * sin(front_angles[i+1]) + elevation_deg,
            center_lat + radius_deg * sin(front_angles[i]) + elevation_deg),
      component = "side",
      station = center_lon,
      segment = i
    )
    side_front <- bind_rows(side_front, quad)
  }
  
  return(list(shadow = shadow, top = top, side = side_front))
}

#################################
# FUNCTION: EXTRACT STOPOVER RASTER IN PLATFORM AREA
#################################

extract_platform_stopover_raster <- function(stopover_raster, center_lon, center_lat, radius_km, elevation_deg) {
  
  # Create circular buffer around radar
  center_sf <- st_sfc(st_point(c(center_lon, center_lat)), crs = 4326)
  buffer_sf <- st_buffer(center_sf, dist = radius_km * 1000)  # Convert km to m
  
  # Crop and mask stopover raster to this circle
  stopover_crop <- crop(stopover_raster, vect(buffer_sf))
  stopover_masked <- mask(stopover_crop, vect(buffer_sf))
  
  # Convert to dataframe
  stopover_df <- as.data.frame(stopover_masked, xy = TRUE, na.rm = TRUE)
  
  if(nrow(stopover_df) == 0) return(NULL)
  
  names(stopover_df)[3] <- "stopover"
  
  # Elevate Y coordinates to platform height
  stopover_df$y <- stopover_df$y + elevation_deg
  stopover_df$station <- paste(center_lon, center_lat, sep="_")
  
  return(stopover_df)
}

#################################
# FUNCTION: EXTRACT DATA FOR YEAR
#################################

extract_year_data <- function(stopover_stack, spei_stack, year) {
  
  cat("  Extracting data for year", year, "\n")
  
  # Extract year from layer names
  stopover_layer_names <- names(stopover_stack)
  stopover_years <- as.numeric(gsub(".*_(\\d{4})_.*", "\\1", stopover_layer_names))
  stopover_idx <- which(stopover_years == year)
  
  spei_layer_names <- names(spei_stack)
  spei_years <- as.numeric(gsub(".*_(\\d{4})_.*", "\\1", spei_layer_names))
  spei_idx <- which(spei_years == year)
  
  if(length(stopover_idx) == 0 || length(spei_idx) == 0) {
    stop("Year ", year, " not found in data")
  }
  
  cat("    Using stopover layer", stopover_idx, "\n")
  cat("    Using SPEI layer", spei_idx, "\n")
  
  stopover_layer <- stopover_stack[[stopover_idx]]
  spei_layer <- spei_stack[[spei_idx]]
  
  stopover_layer <- project(stopover_layer, "EPSG:4326")
  spei_layer <- project(spei_layer, "EPSG:4326")
  
  return(list(stopover = stopover_layer, spei = spei_layer))
}

#################################
# FUNCTION: CREATE PANEL PLOT
#################################

create_panel <- function(stopover_raster, spei_raster, panel_title, panel_label) {
  
  cat("  Creating elevated platforms for", nrow(radar_stations), "radars...\n")
  
  # Build platform geometry for all radars
  all_shadows <- data.frame()
  all_sides <- data.frame()
  all_tops <- data.frame()
  
  for(i in 1:nrow(radar_stations)) {
    platform <- create_platform_geometry(
      center_lon = radar_stations$lon[i],
      center_lat = radar_stations$lat[i],
      radius_km = platform_radius_km,
      elevation_deg = platform_elevation_deg
    )
    
    all_shadows <- bind_rows(all_shadows, platform$shadow)
    all_sides <- bind_rows(all_sides, platform$side)
    all_tops <- bind_rows(all_tops, platform$top)
  }
  
  # Extract stopover data on platforms as raster tiles
  cat("  Extracting stopover patterns on platforms...\n")
  all_stopover_raster <- data.frame()
  
  for(i in 1:nrow(radar_stations)) {
    raster_df <- extract_platform_stopover_raster(
      stopover_raster,
      center_lon = radar_stations$lon[i],
      center_lat = radar_stations$lat[i],
      radius_km = platform_radius_km,
      elevation_deg = platform_elevation_deg
    )
    
    if(!is.null(raster_df)) {
      all_stopover_raster <- bind_rows(all_stopover_raster, raster_df)
    }
  }
  
  # Normalize stopover for visualization
  if(nrow(all_stopover_raster) > 0) {
    all_stopover_raster$stopover_norm <- all_stopover_raster$stopover / max(all_stopover_raster$stopover, na.rm = TRUE)
  }
  
  # Prepare SPEI background
  cat("  Preparing SPEI background...\n")
  tx_la_buffered <- st_buffer(tx_la, dist = 0.5)
  spei_masked <- mask(spei_raster, vect(tx_la_buffered))
  spei_df <- as.data.frame(spei_masked, xy = TRUE, na.rm = TRUE)
  names(spei_df)[3] <- "spei"
  
  # Build plot
  cat("  Building plot...\n")
  
  p <- ggplot() +
    # SPEI drought background FIRST (at bottom)
    geom_raster(data = spei_df, aes(x = x, y = y, fill = spei)) +
    scale_fill_gradientn(
      colors = c("#8B4513", "#CD853F", "#DEB887", "#F5DEB3", "#F0F8FF", "#B0E0E6", "#4682B4", "#1E90FF"),
      name = "SPEI\n(Drought Index)",
      limits = c(-2.5, 2.5),
      na.value = NA,
      breaks = c(-2, -1, 0, 1, 2),
      labels = c("-2 Extreme Drought", "-1 Dry", "0 Normal", "1 Wet", "2 Very Wet")
    ) +
    # State boundaries
    geom_sf(data = tx_la, fill = NA, color = "gray20", linewidth = 0.8) +
    # BCR boundaries
    geom_sf(data = bcr_clipped, fill = NA, color = "gray50", 
            linewidth = 0.25, linetype = "dotted", alpha = 0.6) +
    # Platform shadows
    geom_polygon(data = all_shadows,
                 aes(x = x, y = y, group = station),
                 fill = "black", color = NA,
                 alpha = platform_shadow_alpha) +
    # Platform sides (cylindrical walls)
    geom_polygon(data = all_sides,
                 aes(x = x, y = y, group = interaction(station, segment)),
                 fill = platform_side_color, 
                 color = "#1a1a1a",
                 alpha = 0.9, linewidth = 0.2) +
    # Stopover patterns ON platforms (as raster tiles - will override SPEI fill)
    {if(nrow(all_stopover_raster) > 0) {
      list(
        ggnewscale::new_scale_fill(),
        geom_raster(data = all_stopover_raster, aes(x = x, y = y, fill = stopover)),
        scale_fill_gradientn(
          colors = c("#fde724", "#7ad151", "#22a884", "#2a788e", "#414487", "#440154"),
          name = "Stopover\nDensity",
          na.value = "transparent"
        )
      )
    }} +
    # North arrow and scale
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
      title = panel_title,
      subtitle = panel_label,
      x = "Longitude",
      y = "Latitude"
    ) +
    coord_sf(xlim = c(-108.5, -88.5), ylim = c(24, 37), 
             crs = flat_crs, expand = FALSE) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 11, hjust = 0, color = "gray30"),
      legend.position = "right",
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank()
    )
  
  return(p)
}

#################################
# GENERATE PANELS
#################################

cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("GENERATING PANELS\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# Panel (a): Fall Drought Year
cat("Panel (a): Fall", fall_drought_year, "(Drought)\n")
fall_drought_data <- extract_year_data(fall_stopover, fall_spei, fall_drought_year)
p_fall_drought <- create_panel(
  stopover_raster = fall_drought_data$stopover,
  spei_raster = fall_drought_data$spei,
  panel_title = paste0("(a) Fall ", fall_drought_year),
  panel_label = "Wet Year"
)

print(p_fall_drought)

ggsave(
  filename = file.path(output_dir, paste0("panel_a_fall_", fall_drought_year, "_drought.png")),
  plot = p_fall_drought,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi
)
cat("  ✓ Panel (a) saved\n\n")

# Add similar blocks for panels b, c, d...

# Panel (b): Spring Drought Year

cat("Panel (b): Spring", spring_drought_year, "(Drought)\n")
spring_drought_data <- extract_year_data(spring_stopover, spring_spei, spring_drought_year)
p_spring_drought <- create_panel(
  stopover_raster = spring_drought_data$stopover,
  spei_raster = spring_drought_data$spei,
  panel_title = paste0("(b) Spring ", spring_drought_year),
  panel_label = "Drought Year"
)

print(p_spring_drought)

ggsave(
  filename = file.path(output_dir, paste0("panel_b_spring_", spring_drought_year, "_drought.png")),
  plot = p_spring_drought,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi
)

cat(paste(rep("=", 70), collapse=""), "\n")
cat("FIGURE 3 COMPLETE!\n")
cat(paste(rep("=", 70), collapse=""), "\n")