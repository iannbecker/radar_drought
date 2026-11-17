#################################
#
#   Figure 5 - SPEI Main Effect, Fall and Spring
#   Ian Becker
#   November 2026
#
################################

# Partial effect of the winning SPEI window for each season,
# trimmed to the 5th-95th percentile of observed SPEI, on a shared y-axis

library(mgcv)
library(gratia)
library(ggplot2)
library(dplyr)

output_dir <- "figures_tables"
dir.create(output_dir, showWarnings = FALSE)

# EXTRACT SMOOTHS ------------------------------------

get_spei <- function(mod, term, var) {
  sm <- smooth_estimates(mod, select = term, n = 200) %>% add_confint()
  q  <- quantile(mod$model[[var]], c(0.05, 0.95))
  sm <- sm[sm[[var]] >= q[1] & sm[[var]] <= q[2], ]
  sm$spei <- sm[[var]]
  sm
}

# Fall: current-season window
fall_mod <- readRDS("fall_BEST_FINAL_full_interaction.rds")
fall_sm  <- get_spei(fall_mod, "s(SPEI_fall_current)", "SPEI_fall_current")
rm(fall_mod); gc()

# Spring: previous-year same-season window
spring_mod <- readRDS("spring_BEST_FINAL_full_interaction.rds")
spring_sm  <- get_spei(spring_mod, "s(SPEI_spring_prev_year)", "SPEI_spring_prev_year")
rm(spring_mod); gc()

# If SPEI was scaled before modeling, back-transform here with saved scaling info:
# fall_sm$spei   <- fall_sm$spei   * fall_sd   + fall_mean
# spring_sm$spei <- spring_sm$spei * spring_sd + spring_mean

# SHARED Y-AXIS ------------------------------------

y_lim <- range(c(fall_sm$.lower_ci, fall_sm$.upper_ci,
                 spring_sm$.lower_ci, spring_sm$.upper_ci))

# PLOT FUNCTION ------------------------------------

plot_spei <- function(sm, x_lab) {
  ggplot(sm, aes(x = spei, y = .estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), fill = "#D55E00", alpha = 0.5) +
    geom_line(linewidth = 1, color = "#D55E00") +
    coord_cartesian(ylim = y_lim) +
    labs(x = x_lab, y = "Partial effect (log scale)") +
    theme_classic(base_size = 15)
}

p_fall   <- plot_spei(fall_sm,   "SPEI (current season)")
p_spring <- plot_spei(spring_sm, "SPEI (previous year, same season)")

# SAVE ------------------------------------

ggsave(file.path(output_dir, "figure5a_spei_fall.png"),   p_fall,   width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "figure5b_spei_spring.png"), p_spring, width = 6, height = 5, dpi = 300)


make_wheel <- function(mig_months, win_months, win_year, col) {
  base <- expand.grid(month = 1:12, ring = c("previous", "current"))
  base$ymin <- ifelse(base$ring == "previous", 1, 2.2)
  base$ymax <- base$ymin + 1
  base$fill <- ifelse(base$ring == win_year & base$month %in% win_months, "window", "none")
  
  ggplot() +
    geom_rect(data = base,
              aes(xmin = month - 0.5, xmax = month + 0.5, ymin = ymin, ymax = ymax, fill = fill),
              color = "white", linewidth = 0.6) +
    annotate("rect", xmin = min(mig_months) - 0.5, xmax = max(mig_months) + 0.5,
             ymin = 2.2, ymax = 3.2, fill = NA, color = "black", linewidth = 0.9) +  # migration
    annotate("text", x = 1:12, y = 3.75, label = substr(month.abb, 1, 1), size = 3.5) +
    scale_fill_manual(values = c(none = "grey90", window = col), guide = "none") +
    scale_x_continuous(limits = c(0.5, 12.5)) +
    scale_y_continuous(limits = c(0, 4)) +
    coord_polar(start = 0) +   # January at 12 o'clock
    theme_void()
}

# placeholders: fall migration Aug–Nov, window = current season
wheel_fall   <- make_wheel(mig_months = 8:11, win_months = 8:11, win_year = "current",  col = "#D55E00")
# placeholders: spring migration Mar–May, window = previous year same season
wheel_spring <- make_wheel(mig_months = 3:5,  win_months = 3:5,  win_year = "previous", col = "#56B4E9")

ggsave("figures_tables/figure5_wheel_fall.png",   wheel_fall,   width = 2.5, height = 2.5, dpi = 300, bg = "white")
ggsave("figures_tables/figure5_wheel_spring.png", wheel_spring, width = 2.5, height = 2.5, dpi = 300, bg = "white")
