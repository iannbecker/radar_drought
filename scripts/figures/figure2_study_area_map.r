#################################
#
#   Figure 2 - Study Area Map
#   Ian Becker
#   November 2025
#
################################

library(sf)
library(ggplot2)
library(ggspatial)
library(dplyr)
library(tigris)
library(viridis)

options(tigris_use_cache = TRUE)

# BCR lookup table

bcr_lookup <- data.frame(
  code = c(1, 2, 4, 5, 6, 7, 9, 12, 14, 15),
  name = c(
    "Central Mixed Grass Prairie", "Chihuahuan Desert",
    "Edwards Plateau", "Gulf Coastal Prairie",
    "Mississippi Alluvial Valley", "Oaks And Prairies",
    "Shortgrass Prairie", "Southeastern Coastal Plain",
    "Tamaulipan Brushlands", "West Gulf Coastal Plain/Ouachitas"
  ),
  stringsAsFactors = FALSE
)

#################################
# LOAD DATA
#################################

cat("Downloading state boundaries...\n")
texas <- states(cb = TRUE) %>% filter(NAME == "Texas")
louisiana <- states(cb = TRUE) %>% filter(NAME == "Louisiana")
study_area <- bind_rows(texas, louisiana)

# Load BCR shapefile

gdb_path <- "/Users/ianbecker/Library/CloudStorage/OneDrive-TheUniversityofTexas-RioGrandeValley/DroughtRadar/NABCI_ecoregion.gdb"
bcr_shapefile <- st_read(gdb_path, layer = "BCR_Terrestrial_Master")

# Transform to same CRS

bcr_shapefile <- st_transform(bcr_shapefile, st_crs(study_area))

# Crop BCRs to study area

cat("Cropping BCRs to Texas and Louisiana...\n")
sf_use_s2(FALSE)
bcr_shapefile <- st_make_valid(bcr_shapefile)
study_area <- st_make_valid(study_area)
bcr_shapefile <- st_intersection(bcr_shapefile, st_union(study_area))
sf_use_s2(TRUE)

# Find BCR column name and join with lookup table

bcr_id_col <- NULL
for(col in c("bcr", "BCR", "bcr_code", "BCR_CODE", "bcr_label_name", "BCR_LABEL_NAME")) {
  if(col %in% names(bcr_shapefile)) {
    bcr_id_col <- col
    break
  }
}

cat("Using BCR column:", bcr_id_col, "\n")

# Join with lookup table to get clean names

if(is.numeric(bcr_shapefile[[bcr_id_col]])) {
  bcr_shapefile <- bcr_shapefile %>%
    left_join(bcr_lookup, by = setNames("code", bcr_id_col))
} else {
  bcr_shapefile <- bcr_shapefile %>%
    left_join(bcr_lookup, by = setNames("name", bcr_id_col))
}

# Calculate centroids for labels

bcr_centroids <- st_centroid(bcr_shapefile)

#################################
# RADAR DATA AND DETECTION CIRCLES
#################################

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

#################################
# CREATE MAP
#################################

cat("Creating map...\n")

p <- ggplot() +
  geom_sf(data = study_area, fill = NA, color = "black", linewidth = 0.8) + # state boundaries
  geom_sf(data = bcr_shapefile, aes(fill = bcr_label_name), 
          color = "white", linewidth = 0.5, alpha = 0.5) + # ecoregions
  geom_sf(data = radar_detection_range, fill = "red", alpha = 0.1, 
          color = "red", linewidth = 0.5) + # Radar detection radii
  geom_sf_text(data = radar_locations, aes(label = site), 
              size = 3, fontface = "bold") + # Radar labels
  scale_fill_viridis_d(option = "turbo", name = "Ecoregion") +
  coord_sf(datum = st_crs(study_area)) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 15),
    panel.grid = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)
  ) +
  annotation_scale(location = "tr", width_hint = 0.3) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  guides(fill = guide_legend(
    nrow = 5, 
    keyheight = unit(0.4, "cm"),
    keywidth = unit(0.4, "cm")
  ))

# Save figure

output_path <- "figures_tables/figure2_study_area_map.png"
ggsave(output_path, plot = p, width = 10, height = 8, dpi = 300)
