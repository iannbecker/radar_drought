##############################
#
#   Figure 3 Part 1: Stopover density
#   Ian Becker
#   11/12/2025
#
##############################

library(terra)
library(sf)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(ggspatial)  
library(tigris)

options(tigris_use_cache = TRUE)
setwd("PATH HERE")

#################################
# DATA FORMATTING
#################################

# Load stopover rasters

fall_stopover <- rast("PATH HERE")
spring_stopover <- rast("PATH HERE")

# Calculate mean across all years

fall_mean <- mean(fall_stopover, na.rm = TRUE)
spring_mean <- mean(spring_stopover, na.rm = TRUE)

# Load BCR shapefile (for boundaries only)

bcr <- st_read("NABCI_ecoregion.gdb", layer = "BCR_Terrestrial_Master", quiet = TRUE)

# Make sure CRS matches

bcr <- st_transform(bcr, crs(fall_mean))

# Crop BCR to study area

study_extent <- ext(fall_mean)
bcr_crop <- st_crop(bcr, study_extent)
nrow(bcr_crop) # Check number of BCRS

# Get state boundaries using tigris

states <- states(cb = TRUE, resolution = "20m")

# Get just Texas and Louisiana (keep separate, don't union)

tx_la <- states %>%
  filter(STUSPS %in% c("TX", "LA")) %>%
  st_transform(crs(fall_mean))

# Lat/long projection for flat mapping

flat_crs <- st_crs(4326)  # WGS84 lat/lon

# Reproject everything to flat projection

tx_la <- st_transform(tx_la, flat_crs)
bcr_crop <- st_transform(bcr_crop, flat_crs)

# Reproject rasters

fall_mean <- project(fall_mean, "EPSG:4326")
spring_mean <- project(spring_mean, "EPSG:4326")

# Clip BCR to only show within TX/LA boundaries

bcr_clipped <- st_intersection(bcr_crop, tx_la)
plot(bcr_clipped) # Check clipping

cat("Data prep complete!\n")

#################################
# RADAR STATION LOCATIONS
#################################

# Radar station coordinates

texas_radars <- data.frame(
  site = c("KGRK", "KHGX", "KCRP", "KDFX", "KFWS", 
           "KDYX", "KEPZ", "KLBB", "KMAF", "KSJT", "KAMA", "KBRO"),
  name = c("Ft. Hood", "Houston", "Corpus Christi", 
           "Del Rio", "Dallas/Fort Worth", "Dyess AFB", "El Paso", 
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

# Combine all radars

radar_stations <- bind_rows(texas_radars, louisiana_radars) %>%
  rename(station = site, lat = latitude, lon = longitude)

# Convert to sf object - coordinates are in WGS84 (lat/lon)

radar_sf <- st_as_sf(radar_stations, coords = c("lon", "lat"), crs = 4326)  # WGS84

# Transform to flat projection

radar_sf <- st_transform(radar_sf, flat_crs)

# Extract mean stopover values at each radar location

buffer_dist <- 80000  # 80km i

# Buffer in the raster's CRS
# 80km = 80000 meters

radar_buffered <- st_buffer(radar_sf, dist = buffer_dist)

# Extract values

fall_values <- terra::extract(fall_mean, vect(radar_buffered), fun = mean, na.rm = TRUE)
spring_values <- terra::extract(spring_mean, vect(radar_buffered), fun = mean, na.rm = TRUE)

# Check for NAs

cat("\nFall NAs:", sum(is.na(fall_values[[2]])), "of", nrow(fall_values), "\n")
cat("Spring NAs:", sum(is.na(spring_values[[2]])), "of", nrow(spring_values), "\n")

# Add to dataframe with ORIGINAL lat/lon preserved for plotting

radar_data <- radar_stations %>%
  mutate(
    fall_stopover = fall_values[[2]],  
    spring_stopover = spring_values[[2]]
  )

# Convert to sf for easier plotting with transformed coordinates

radar_data_sf <- st_as_sf(radar_data, coords = c("lon", "lat"), crs = 4326)
radar_data_sf <- st_transform(radar_data_sf, flat_crs)

# Add the transformed coordinates back to the dataframe

coords <- st_coordinates(radar_data_sf)
radar_data$x <- coords[,1]
radar_data$y <- coords[,2]

cat("Radar data locations prepped!\n")

#################################
# CREATE MAPS
#################################

# Get study area boundary for plotting

study_bbox <- st_as_sfc(st_bbox(bcr_crop))

# Color palette for stopover (yellow to red)

stopover_colors <- scale_fill_gradientn(
  colors = c("#FFF7BC", "#FEE391", "#FEC44F", "#FE9929", "#EC7014", "#CC4C02", "#8C2D04"),
  name = "Mean\nStopover\nDensity",
  na.value = "gray90"
)

stopover_color_scale <- scale_color_gradientn(
  colors = c("#FFF7BC", "#FEE391", "#FEC44F", "#FE9929", "#EC7014", "#CC4C02", "#8C2D04"),
  name = "Mean\nStopover\nDensity",
  na.value = "gray90"
)

# Set up common range for both panels

value_range <- range(c(radar_data$fall_stopover, radar_data$spring_stopover), na.rm = TRUE)

#################################
# FALL
#################################

p_fall <- ggplot() +
  geom_sf(data = tx_la, fill = "gray95", color = "gray40", linewidth = 0.5) + # State boundaries
  geom_sf(data = bcr_clipped, fill = NA, color = "gray60", linewidth = 0.4) + # BCR boundaries
  geom_sf(data = radar_buffered, aes(fill = radar_data$fall_stopover), # Radar stations
          color = "black", linewidth = 0.5, alpha = 0.8) +
  scale_fill_gradientn(
    colors = c("#FFF7BC", "#FEE391", "#FEC44F", "#FE9929", "#EC7014", "#CC4C02", "#8C2D04"),
    name = "Mean Stopover",
    limits = NULL
  ) +
  annotation_north_arrow(
    location = "tr", 
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1.5, "cm"),
    width = unit(1.5, "cm")
  ) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  labs(
    title = "(a) Fall Migration",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    panel.grid = element_line(color = "gray90", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  coord_sf(crs = flat_crs, datum = flat_crs)

#################################
# SPRING
#################################

p_spring <- ggplot() +
  geom_sf(data = tx_la, fill = "gray95", color = "gray40", linewidth = 0.5) +
  geom_sf(data = bcr_clipped, fill = NA, color = "gray60", linewidth = 0.4) +
  geom_sf(data = radar_buffered, aes(fill = radar_data$spring_stopover), 
          color = "black", linewidth = 0.5, alpha = 0.8) +
  scale_fill_gradientn(
    colors = c("#FFF7BC", "#FEE391", "#FEC44F", "#FE9929", "#EC7014", "#CC4C02", "#8C2D04"),
    name = "Mean Stopover",
    limits = value_range
  ) +
  labs(
    title = "(b) Spring Migration",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    panel.grid = element_line(color = "gray90", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  coord_sf(crs = flat_crs, datum = flat_crs, expand = TRUE)

# Save plot data for future usage

write.csv(radar_data, 
          file.path(output_dir, "radar_stopover_values.csv"),
          row.names = FALSE)


