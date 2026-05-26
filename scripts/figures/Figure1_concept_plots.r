#################################
#
#   Figure 1 - Study Area Concept Figures 
#   Ian Becker
#   March 2026
#
################################

# This script is used to generate the concept plots used in the Figure 1 visualization

library(ggplot2)
library(dplyr)

# Set seed for reproducibility

set.seed(42)

# RESOURCE AVAILABILITY ------------------------------------

# Sim data

abundance <- data.frame(
  condition = c("Non-Drought", "Drought"),
  arthropod_abundance = c(85, 30)  # hypothetical values
)

# Plot

ggplot(abundance, aes(x = condition, y = arthropod_abundance, fill = condition)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Non-Drought" = "#4a9e6b", "Drought" = "#c95f2e")) +
  labs(
    x = NULL,
    y = "Arthropod Abundance",
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25)
  )

# STOPOVER QUALITY ------------------------------------

# Low-quality habitat

drought <- data.frame(
  precipitation = runif(12, 10, 45),
  mass_gain_rate = runif(12, 5, 30),
  condition = "Drought"
)

# High-quality habitat

non_drought <- data.frame(
  precipitation = runif(12, 35, 80),
  mass_gain_rate = runif(12, 25, 65),
  condition = "Non-Drought"
)

# Combine 

dat <- rbind(drought, non_drought)
dat$condition <- factor(dat$condition, levels = c("Non-Drought", "Drought"))

# Plot

ggplot(dat, aes(x = precipitation, y = mass_gain_rate, color = condition)) +
  geom_point(size = 7, shape = 18) + 
  scale_color_manual(values = c("Non-Drought" = "#4a8c3f", "Drought" = "#c97a2e")) +
  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, 50)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    x = "Precipitation (mm)",
    y = "Daily Mass Gain Rate (%)",
    color = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.75, 0.2),
    legend.text = element_text(size = 20),
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25),
  )

# MASS MORTALITY EVENT PLOT ------------------------------------

# Set parameters

n <- 100
time <- 1:n

# Set trend for event

trend <- ifelse(time <= 65,
                100 - 0.6 * time,
                100 - 0.6 * 65 - 6.5 * (time - 65)
)
trend <- pmax(trend, 12)

# Add noise after event 

trend <- trend + rnorm(n, 0, 4)
trend <- pmax(trend, 2)  # nothing below 0

# Create data frame

dat <- data.frame(time = time, abundance = trend)

# Plot

ggplot(dat, aes(x = time, y = abundance)) +
  geom_line(color = "#cc2222", linewidth = 1.1) +
  scale_x_continuous(
    limits = c(0, n),
    expand = c(0, 0),
    breaks = NULL,
    name = "Time"
  ) +
  scale_y_continuous(
    limits = c(0, 115),
    expand = c(0, 0),
    breaks = NULL,
    name = "Abundance"
  ) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 25, hjust = 0.5, color = "gray30"),
    axis.title.y = element_text(size = 25, color = "gray30"),
    plot.margin = margin(20, 20, 10, 10)
  )

# THERMAL TOLERANCE PLOT ------------------------------------

# Create data frame

temp <- seq(0, 100, length.out = 500)
performance <- 100 * exp(-((temp - 45)^2) / (2 * 15^2))

dat <- data.frame(temp = temp, performance = performance)

# Shading regions

dat$region <- case_when(
  dat$temp < 20                        ~ "cold",
  dat$temp >= 20 & dat$temp <= 70      ~ "optimal",
  dat$temp > 70                        ~ "exceeded"
)

drought_line <- 75  # where drought temp sits 

# Plot

ggplot(dat, aes(x = temp, y = performance)) +
  
  # Shaded regions
  
  geom_area(data = subset(dat, region == "exceeded"),
            fill = "#e83030", alpha = 0.2) +
  geom_area(data = subset(dat, region == "optimal"),
            fill = "#4a8c3f", alpha = 0.12) +
  
  # Main curve
  
  geom_line(linewidth = 1.3, color = "black") +
  
  # Drought threshold line
  
  geom_vline(xintercept = drought_line, linetype = "dashed",
             color = "#e83030", linewidth = 0.9) +
  
  # Labels
  
  annotate("text", x = drought_line + 2, y = 85,
           label = "Drought\nTemp", color = "#e83030",
           hjust = 0, size = 6, fontface = "bold") +
  annotate("text", x = 45, y = 106,
           label = "Optimal", color = "#4a8c3f",
           hjust = 0.5, size = 6) +
  scale_x_continuous(
    breaks = c(10, 45, 75),
    labels = c("Cold", "Optimal Temp", "Critical\nThreshold"),
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(limits = c(0, 115), expand = c(0, 0)) +
  labs(x = "Temperature", y = "Performance") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 15, color = "gray20"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size = 25),
    plot.margin = margin(15, 15, 10, 10)
  )
