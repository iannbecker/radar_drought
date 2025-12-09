##############################
#
#   Drought timeseries - Figure 2
#   Ian Becker
#   11/20/2025
#
##############################

library(terra)
library(ggplot2)
library(dplyr)
library(sf)
library(tidyr)
library(tigris)

options(tigris_use_cache = TRUE)
setwd("PATH HERE")

#################################
# LOAD AND CLEAN SPEI DATA
#################################

############### Prep SPEI data

# Note this is starting with filtering the entire dataset - should be saved hereafter

# Load SPEI full dataset

spei <- rast("nclimgrid-spei-pearson-01.nc")

# Filter to 1995-2020

time <- time(spei)
years<- 1995:2020
time_indices <- which(format(time, "%Y") %in% years)
spei_1995_2020 <- spei[[time_indices]]

# Save raster for future use

writeRaster(spei_1995_2020, "spei_FULL_1995_2020.tif", overwrite=TRUE)

############### Prep study area

# Load in raster data

spei <- rast("spei_FULL_1995_2020.tif")

# Load TX and LA state boundaries

tx_la <- states(cb = TRUE) %>%
  filter(STUSPS %in% c("TX", "LA"))

# Transform to match SPEI CRS

tx_la <- st_transform(tx_la, crs(spei))

# Crop SPEI to study area

cat("Cropping SPEI to TX/LA study area...\n")
spei_study <- crop(spei, tx_la)
spei_study <- mask(spei_study, tx_la)

#################################
# CATEGORIZE DROUGHT
#################################

cat("Categorizing drought severity for each month...\n")

# Define drought categories (standard SPEI thresholds)

categorize_drought <- function(spei_values) {
  categories <- case_when(
    spei_values <= -2.0 ~ "Extremely Dry",
    spei_values <= -1.5 ~ "Very Dry",
    spei_values <= -1.0 ~ "Moderately Dry",
    spei_values <= -0.5 ~ "Abnormally Dry",
    TRUE ~ "No Drought"
  )
  return(categories)
}

# Get time stamps

dates <- time(spei_study)

# Initialize results dataframe

drought_history <- data.frame()

# Process each month

for(i in 1:nlyr(spei_study)) {
  if(i %% 50 == 0) cat("  Processing layer", i, "of", nlyr(spei_study), "\n")
  
  # Get SPEI values for this month
  spei_vals <- values(spei_study[[i]])
  spei_vals <- spei_vals[!is.na(spei_vals)]  # Remove NAs
  
  # Categorize
  categories <- categorize_drought(spei_vals)
  
  # Calculate percentages
  total_pixels <- length(spei_vals)
  pct_table <- table(categories) / total_pixels * 100
  
  # Store results
  month_data <- data.frame(
    date = dates[i],
    category = names(pct_table),
    percent = as.numeric(pct_table),
    stringsAsFactors = FALSE
  )
  
  drought_history <- bind_rows(drought_history, month_data)
}

cat("\nDrought categorization complete!\n")
cat("Total observations:", nrow(drought_history), "\n\n")

#################################
# PREPARE FOR PLOTTING
#################################

cat("Preparing data for plotting...\n")

# Set factor levels (from most severe to least)

drought_history$category <- factor(
  drought_history$category,
  levels = c("Extremely Dry", "Very Dry", "Moderately Dry", 
             "Abnormally Dry", "No Drought")
)

# Extract year and season for potential down the road filtering

drought_history$year <- as.numeric(format(drought_history$date, "%Y"))
drought_history$month <- as.numeric(format(drought_history$date, "%m"))

# Add season labels for potential down the road use

drought_history$season <- case_when(
  drought_history$month %in% c(9, 10) ~ "Fall",
  drought_history$month %in% c(3, 4, 5) ~ "Spring",
  TRUE ~ "Other"
)

#################################
# CREATE STACKED AREA PLOT
#################################

cat("Creating stacked area plot...\n")

# Define colors 

drought_colors <- c(
  "Extremely Dry" = "#730000",  # Dark red/maroon
  "Very Dry" = "#E60000",      # Bright red
  "Moderately Dry" = "#FFA500",       # Orange
  "Abnormally Dry" = "#FFFF00", # Yellow
  "No Drought" = "white"           # White/none
)

drought_history_plot <- drought_history %>%
  filter(category != "No Drought")

p <- ggplot(drought_history_plot, aes(x = date, y = percent, fill = category)) +
  geom_area(position = "stack", alpha = 0.9) +
  scale_fill_manual(
    values = drought_colors,
    name = NULL,
    breaks = c("Extremely Dry", "Very Dry", "Moderately Dry", "Abnormally Dry")
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0),
    breaks = seq(0, 100, 20)
  ) +
  scale_x_continuous(
    breaks = seq(as.Date("1995-01-01"), as.Date("2020-12-01"), by = "5 years"),
    labels = seq(1995, 2020, 5)
  ) +
  labs(
    x = "Year",
    y = "Percent of Study Area in Drought"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(nrow = 1, reverse = TRUE))

# Save plot

output_dir <- "gam_results/drought_figures"
output_file <- file.path(output_dir, "drought_history_stacked.png")
ggsave(output_file, p, width = 12, height = 6, dpi = 300)
cat("\nPlot saved:", output_file, "\n")

cat("All done!\n")