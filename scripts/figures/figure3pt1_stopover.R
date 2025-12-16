##############################
#
#   Figure 3: 3D Stopover Bars on Drought Maps
#   Ian Becker
#   12/16/2025
#
#   Creates individual panels showing stopover density as 3D bars
#   overlaid on drought severity maps for exemplar years
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
library(scales)

options(tigris_use_cache = TRUE)
setwd("/Users/ianbecker/Library/CloudStorage/OneDrive-TheUniversityofTexas-RioGrandeValley/DroughtRadar")

#################################
# CONFIGURATION - EASY TO CHANGE
#################################

# SELECT EXEMPLAR YEARS
fall_drought_year <- 2011  # Extremely dry fall
fall_wet_year <- 2015      # Wet fall
spring_drought_year <- 2006  # Dry spring
spring_wet_year <- 2010     # Wet spring

# File paths
fall_stopover_file <- "raster_data/stopover_fall.tif"
spring_stopover_file <- "raster_data/stopover_spring.tif"
fall_spei_file <- "raster_data/cropped_rasters/fall_raster/SPEI_1995_2020.tif"  
spring_spei_file <- "raster_data/cropped_rasters/spring_raster/SPEI_spring_1995_2020.tif"
bcr_gdb <- "NABCI_ecoregion.gdb"

# Output settings
output_dir <- "figure3_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 3D bar plot settings
bar_width_degrees <- 0.7    # Width of bars in decimal degrees (~70 km) - INCREASED
bar_height_scale <- 1.5     # Multiplier for bar height - INCREASED
bar_3d_offset <- 0.25       # 3D perspective offset in degrees
bar_color <- "#34495e"      # Dark slate blue-gray (main front face)
bar_color_top <- "#4a6b87"  # Lighter for top face (catching light)
bar_color_right <- "#1a252f" # Darkest for right face (shadow side)
bar_edge_color <- "#0d1419" # Very dark edge for contrast
bar_alpha <- 0.9
bar_shadow_alpha <- 0.3     # Shadow transparency for depth effect

# Plot dimensions
plot_width <- 8
plot_height <- 7
plot_dpi <- 300

cat(paste(rep("=", 70), collapse=""), "\n")
cat("FIGURE 3: 3D STOPOVER BARS ON DROUGHT MAPS\n")
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

cat("  Stopover data: ", nlyr(fall_stopover), "years (", min(years), "-", max(years), ")\n", sep="")
cat("  SPEI data: ", nlyr(fall_spei), "layers\n", sep="")

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

radar_sf <- st_as_sf(radar_stations, coords = c("lon", "lat"), crs = 4326)
radar_stations$lon_orig <- radar_stations$lon
radar_stations$lat_orig <- radar_stations$lat

cat("  Radar stations:", nrow(radar_stations), "\n\n")

#################################
# FUNCTION: CREATE 3D BAR POLYGONS
#################################

create_3d_bar <- function(x, y, height, width = bar_width_degrees, 
                          offset = bar_3d_offset) {
  h <- height * bar_height_scale
  
  # Shadow (on ground, behind bar)
  shadow_offset <- 0.08
  shadow <- data.frame(
    x = c(x - width/2, x + width/2, x + width/2 + offset*1.5, x - width/2 + offset*1.5, x - width/2),
    y = c(y - shadow_offset, y - shadow_offset, y - shadow_offset + offset*0.5, 
          y - shadow_offset + offset*0.5, y - shadow_offset),
    face = "shadow"
  )
  
  # Front face
  front <- data.frame(
    x = c(x - width/2, x + width/2, x + width/2, x - width/2, x - width/2),
    y = c(y, y, y + h, y + h, y),
    face = "front"
  )
  
  # Top face (parallelogram)
  top <- data.frame(
    x = c(x - width/2, x + width/2, x + width/2 + offset, x - width/2 + offset, x - width/2),
    y = c(y + h, y + h, y + h + offset, y + h + offset, y + h),
    face = "top"
  )
  
  # Right face (parallelogram)
  right <- data.frame(
    x = c(x + width/2, x + width/2 + offset, x + width/2 + offset, x + width/2, x + width/2),
    y = c(y, y + offset, y + h + offset, y + h, y),
    face = "right"
  )
  
  bar_data <- bind_rows(shadow, front, top, right)
  bar_data$x_center <- x
  bar_data$y_center <- y
  bar_data$height <- height
  
  return(bar_data)
}

#################################
# FUNCTION: EXTRACT DATA FOR SPECIFIC YEAR
#################################


extract_year_data <- function(stopover_stack, spei_stack, year, years_vector) {
  
  cat("  Extracting data for year", year, "\n")
  
  # Extract year from layer names instead of using positional index
  # Layer names are like "VIR_2011_fall" or "SPEI_2011_fall"
  
  # For stopover stack
  stopover_layer_names <- names(stopover_stack)
  stopover_years <- as.numeric(gsub(".*_(\\d{4})_.*", "\\1", stopover_layer_names))
  stopover_idx <- which(stopover_years == year)
  
  # For SPEI stack
  spei_layer_names <- names(spei_stack)
  spei_years <- as.numeric(gsub(".*_(\\d{4})_.*", "\\1", spei_layer_names))
  spei_idx <- which(spei_years == year)
  
  if(length(stopover_idx) == 0) {
    stop("Year ", year, " not found in stopover data\n",
         "Available years: ", paste(range(stopover_years), collapse="-"))
  }
  
  if(length(spei_idx) == 0) {
    stop("Year ", year, " not found in SPEI data\n",
         "Available years: ", paste(range(spei_years), collapse="-"))
  }
  
  cat("    Stopover: layer", stopover_idx, "(", stopover_layer_names[stopover_idx], ")\n")
  cat("    SPEI: layer", spei_idx, "(", spei_layer_names[spei_idx], ")\n")
  
  stopover_layer <- stopover_stack[[stopover_idx]]
  spei_layer <- spei_stack[[spei_idx]]
  
  stopover_layer <- project(stopover_layer, "EPSG:4326")
  spei_layer <- project(spei_layer, "EPSG:4326")
  
  buffer_dist <- 80000  
  
  radar_sf_original <- st_as_sf(radar_stations, 
                                coords = c("lon_orig", "lat_orig"), 
                                crs = 4326)
  radar_buffered <- st_buffer(radar_sf_original, dist = buffer_dist)
  
  stopover_values <- terra::extract(stopover_layer, 
                                    vect(radar_buffered), 
                                    fun = mean, 
                                    na.rm = TRUE)
  
  radar_data <- radar_stations %>%
    mutate(
      stopover = stopover_values[[2]],
      lon = lon_orig,
      lat = lat_orig
    ) %>%
    filter(!is.na(stopover)) %>%
    mutate(stopover_norm = stopover / max(stopover, na.rm = TRUE))
  
  # DIAGNOSTIC: Check if Brownsville is included
  brownsville_included <- "KBRO" %in% radar_data$station
  cat("  Brownsville (KBRO) included:", brownsville_included, "\n")
  if(!brownsville_included) {
    cat("  WARNING: Brownsville radar was filtered out (NA stopover value)\n")
    cat("  This might mean the stopover raster doesn't extend far enough south\n")
  } else {
    bro_data <- radar_data[radar_data$station == "KBRO", ]
    cat("  Brownsville stopover value:", round(bro_data$stopover, 2), 
        "(normalized:", round(bro_data$stopover_norm, 2), ")\n")
  }
  
  return(list(
    stopover = stopover_layer,
    spei = spei_layer,
    radar_data = radar_data
  ))
}

#################################
# FUNCTION: CREATE SINGLE PANEL PLOT
#################################

create_panel <- function(spei_raster, radar_data, panel_title, panel_label) {
  
  cat("  Creating 3D bars for", nrow(radar_data), "radars...\n")
  
  all_bars <- data.frame()
  
  for(i in 1:nrow(radar_data)) {
    bar <- create_3d_bar(
      x = radar_data$lon[i],
      y = radar_data$lat[i],
      height = radar_data$stopover_norm[i]
    )
    bar$station <- radar_data$station[i]
    all_bars <- bind_rows(all_bars, bar)
  }
  
  cat("  Converting SPEI raster to dataframe and masking to state boundaries...\n")
  
  # Buffer state boundaries to ensure we don't clip radar locations near edges
  # This is especially important for Brownsville which is near the southern tip
  tx_la_buffered <- st_buffer(tx_la, dist = 0.5)  # 0.5 degree buffer (~55 km)
  
  # Mask SPEI to buffered state boundaries
  spei_masked <- mask(spei_raster, vect(tx_la_buffered))
  spei_df <- as.data.frame(spei_masked, xy = TRUE, na.rm = TRUE)
  names(spei_df)[3] <- "spei"
  
  cat("  Building plot with enhanced 3D effects...\n")
  
  p <- ggplot() +
    # Background SPEI drought map (now masked to states)
    geom_raster(data = spei_df, aes(x = x, y = y, fill = spei)) +
    scale_fill_gradientn(
      colors = c("#8B4513", "#CD853F", "#DEB887", "#F5DEB3", "#F0F8FF", "#B0E0E6", "#4682B4", "#1E90FF"),
      name = "SPEI\n(Drought\nIndex)",
      limits = c(-2.5, 2.5),
      na.value = NA,  # Changed from gray90 to NA so only states show
      breaks = c(-2, -1, 0, 1, 2),
      labels = c("-2\nExtreme\nDrought", "-1\nDry", "0\nNormal", "1\nWet", "2\nVery\nWet")
    ) +
    # State boundaries (thicker for emphasis)
    geom_sf(data = tx_la, fill = NA, color = "gray20", linewidth = 0.8) +
    # BCR boundaries
    geom_sf(data = bcr_clipped, fill = NA, color = "gray50", 
            linewidth = 0.25, linetype = "dotted", alpha = 0.6) +
    # Bar shadows (for depth)
    geom_polygon(data = all_bars %>% filter(face == "shadow"),
                 aes(x = x, y = y, group = interaction(station, face)),
                 fill = "black", color = NA,
                 alpha = bar_shadow_alpha) +
    # 3D bars - right face (darkest - shadow side)
    geom_polygon(data = all_bars %>% filter(face == "right"),
                 aes(x = x, y = y, group = interaction(station, face)),
                 fill = bar_color_right, 
                 color = bar_edge_color,
                 alpha = bar_alpha, linewidth = 0.3) +
    # 3D bars - top face (lightest - catching light)
    geom_polygon(data = all_bars %>% filter(face == "top"),
                 aes(x = x, y = y, group = interaction(station, face)),
                 fill = bar_color_top, 
                 color = bar_edge_color,
                 alpha = bar_alpha, linewidth = 0.3) +
    # 3D bars - front face (main color)
    geom_polygon(data = all_bars %>% filter(face == "front"),
                 aes(x = x, y = y, group = interaction(station, face)),
                 fill = bar_color, 
                 color = bar_edge_color,
                 alpha = bar_alpha, linewidth = 0.5) +
    # Radar location markers (center of circles)
    # REMOVED - no markers at base of bars
    # North arrow
    annotation_north_arrow(
      location = "tr", 
      which_north = "true",
      style = north_arrow_fancy_orienteering,
      height = unit(1.2, "cm"),
      width = unit(1.2, "cm")
    ) +
    # Scale bar
    annotation_scale(location = "bl", width_hint = 0.25, 
                     text_cex = 0.8) +
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
      legend.title = element_text(face = "bold", size = 9, lineheight = 1.1),
      legend.text = element_text(size = 8, lineheight = 0.9),
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.6, "cm"),
      panel.grid = element_line(color = "gray80", linewidth = 0.2),
      panel.background = element_rect(fill = "#f8f9fa", color = NA),  # Very light gray for depth
      plot.background = element_rect(fill = "white", color = NA),
      axis.text = element_text(size = 9),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(p)
}

#################################
# GENERATE ALL PANELS
#################################

cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("GENERATING PANELS\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# Panel (a): Fall Drought Year
cat("Panel (a): Fall", fall_drought_year, "(Drought)\n")
fall_drought_data <- extract_year_data(fall_stopover, fall_spei, 
                                       fall_drought_year, years)
p_fall_drought <- create_panel(
  spei_raster = fall_drought_data$spei,
  radar_data = fall_drought_data$radar_data,
  panel_title = paste0("(a) Fall ", fall_drought_year),
  panel_label = "Drought Year"
)

cat("  Saving panel (a)...\n")
ggsave(
  filename = file.path(output_dir, paste0("panel_a_fall_", fall_drought_year, "_drought.png")),
  plot = p_fall_drought,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi
)
cat("  ✓ Panel (a) saved\n\n")

# Panel (b): Fall Wet Year
cat("Panel (b): Fall", fall_wet_year, "(Wet)\n")
fall_wet_data <- extract_year_data(fall_stopover, fall_spei, 
                                   fall_wet_year, years)
p_fall_wet <- create_panel(
  spei_raster = fall_wet_data$spei,
  radar_data = fall_wet_data$radar_data,
  panel_title = paste0("(b) Fall ", fall_wet_year),
  panel_label = "Wet Year"
)

cat("  Saving panel (b)...\n")
ggsave(
  filename = file.path(output_dir, paste0("panel_b_fall_", fall_wet_year, "_wet.png")),
  plot = p_fall_wet,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi
)
cat("  ✓ Panel (b) saved\n\n")

# Panel (c): Spring Drought Year
cat("Panel (c): Spring", spring_drought_year, "(Drought)\n")
spring_drought_data <- extract_year_data(spring_stopover, spring_spei, 
                                         spring_drought_year, years)
p_spring_drought <- create_panel(
  spei_raster = spring_drought_data$spei,
  radar_data = spring_drought_data$radar_data,
  panel_title = paste0("(c) Spring ", spring_drought_year),
  panel_label = "Drought Year"
)

cat("  Saving panel (c)...\n")
ggsave(
  filename = file.path(output_dir, paste0("panel_c_spring_", spring_drought_year, "_drought.png")),
  plot = p_spring_drought,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi
)
cat("  ✓ Panel (c) saved\n\n")

# Panel (d): Spring Wet Year
cat("Panel (d): Spring", spring_wet_year, "(Wet)\n")
spring_wet_data <- extract_year_data(spring_stopover, spring_spei, 
                                     spring_wet_year, years)
p_spring_wet <- create_panel(
  spei_raster = spring_wet_data$spei,
  radar_data = spring_wet_data$radar_data,
  panel_title = paste0("(d) Spring ", spring_wet_year),
  panel_label = "Wet Year"
)

cat("  Saving panel (d)...\n")
ggsave(
  filename = file.path(output_dir, paste0("panel_d_spring_", spring_wet_year, "_wet.png")),
  plot = p_spring_wet,
  width = plot_width,
  height = plot_height,
  dpi = plot_dpi
)
cat("  ✓ Panel (d) saved\n\n")

#################################
# SAVE SUMMARY DATA
#################################

cat("Saving summary data...\n")

all_radar_data <- bind_rows(
  fall_drought_data$radar_data %>% mutate(season = "Fall", year = fall_drought_year, type = "Drought"),
  fall_wet_data$radar_data %>% mutate(season = "Fall", year = fall_wet_year, type = "Wet"),
  spring_drought_data$radar_data %>% mutate(season = "Spring", year = spring_drought_year, type = "Drought"),
  spring_wet_data$radar_data %>% mutate(season = "Spring", year = spring_wet_year, type = "Wet")
)

write.csv(all_radar_data, 
          file.path(output_dir, "radar_stopover_by_year.csv"),
          row.names = FALSE)

cat("  ✓ Summary data saved\n\n")

#################################
# FINAL SUMMARY
#################################

cat(paste(rep("=", 70), collapse=""), "\n")
cat("FIGURE 3 COMPLETE!\n")
cat(paste(rep("=", 70), collapse=""), "\n")
cat("Output directory:", output_dir, "\n\n")
cat("Files created:\n")
cat("  - panel_a_fall_", fall_drought_year, "_drought.png\n", sep="")
cat("  - panel_b_fall_", fall_wet_year, "_wet.png\n", sep="")
cat("  - panel_c_spring_", spring_drought_year, "_drought.png\n", sep="")
cat("  - panel_d_spring_", spring_wet_year, "_wet.png\n", sep="")
cat("  - radar_stopover_by_year.csv (data summary)\n\n")
cat("Dimensions:", plot_width, "x", plot_height, "inches at", plot_dpi, "DPI\n")
cat("\nReady for assembly in PowerPoint or graphics software!\n")
cat(paste(rep("=", 70), collapse=""), "\n")
