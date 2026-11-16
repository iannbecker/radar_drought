#################################
#
#   Figure 2 - Study Area Map
#   Ian Becker
#   November 2025
#
################################

# This script is used to make Figure 2 in the main text
# which is a map of the study area showing bird conservation regions (BCR) and radar locations

library(sf)
library(ggplot2)
library(ggspatial)
library(dplyr)
library(tigris)
library(viridis)

options(tigris_use_cache = TRUE)

# LOAD DATA ------------------------------------

cat("Downloading state boundaries...\n")
texas <- states(cb = TRUE) %>% filter(NAME == "Texas")
louisiana <- states(cb = TRUE) %>% filter(NAME == "Louisiana")
study_area <- bind_rows(texas, louisiana)

# RADAR LOCATIONS AND DETECTION RADII ------------------------------------

# Texas radars

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

# Louisiana radars

louisiana_radars <- data.frame(
  site = c("KLCH", "KLIX", "KSHV"),
  name = c("Lake Charles", "New Orleans", "Shreveport"),
  latitude = c(30.125, 30.337, 32.451),
  longitude = c(-93.216, -89.826, -93.841)
)


# Combine

all_radars <- rbind(texas_radars, louisiana_radars)

cat("Total radars loaded:", nrow(all_radars), "\n")

# Convert to sf object

radar_locations <- st_as_sf(all_radars, 
                            coords = c("longitude", "latitude"),
                            crs = 4326)  # WGS84

# Transform to match BCR map CRS

radar_locations <- st_transform(radar_locations, st_crs(study_area))

# NEXRAD detection range: 80 km

radar_detection_range <- st_buffer(radar_locations, dist = 80000)   # 80 km in meters

# CREATE STUDY AREA MAP ------------------------------------

cat("Creating map...\n")

p <- ggplot() +
  geom_sf(data = study_area, fill = NA, color = "black", linewidth = 0.8) + 
  geom_sf(data = radar_detection_range, fill = "red", alpha = 0.2,
          color = "red", linewidth = 0.6) + # 80 km detection radii
  geom_sf_text(data = radar_locations, aes(label = site),
               size = 3, fontface = "bold", color = "black") + # radar abbreviations
  coord_sf(datum = st_crs(study_area)) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)
  ) +
  annotation_scale(location = "tr", width_hint = 0.3) +
  annotation_north_arrow(
    location = "tl",
    which_north = "true",
    pad_x = unit(0.3, "cm"),
    pad_y = unit(0.3, "cm"),
    style = north_arrow_fancy_orienteering()
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  )

p

# Save figure

output_path <- "figures_tables/figure2_study_area_map.png"
ggsave(output_path, plot = p, width = 10, height = 8, dpi = 300)
