# ==============================================================================
# Modelling NAAMP Peeper Calling Onset with Lovett (2013) TS3 - Predicted DFC
# ==============================================================================

# Load Packages ----------------------------------------------------------------
library(tidyverse)
library(lme4)

# Data import ------------------------------------------------------------------
## NAAMP Observations
observed <- read_csv("data_processed/naamp_peepers_onset.csv")

## Predicted DFC from ERA5 Land (~9km resolution, ~81km2 )
predicted1 <- read_csv("data_processed/naamp_route_predicted_dfc_1997_2005.csv")
predicted2 <- read_csv("data_processed/naamp_route_predicted_dfc_2006_2015.csv")

## Predicted DFC from Daymet (~1km resolution, ~1km2)
predicted_daymet1 <- read_csv("data_processed/naamp_route_predicted_dfc_daymet_1997_2005.csv")
predicted_daymet2 <- read_csv("data_processed/naamp_route_predicted_dfc_daymet_2006_2015.csv")

# Data prep --------------------------------------------------------------------
## Predicted ERA5 Land
predicted <- bind_rows(predicted1, predicted2) %>%
  distinct(RouteNumber, year, .keep_all = TRUE)

obs_pred_dfc <- observed %>%
  left_join(predicted, join_by(
    RouteNumber == RouteNumber,
    SurveyYear == year
  )) %>%
  select(-intensity) %>%
  janitor::clean_names() %>%
  rename(
    observed_dfc = "onset_doy",
    predicted_dfc = "first"
  ) %>%
  mutate(
    # Create ordered latitude bins safely
    lat_bin = cut(
      lat,
      breaks = c(-Inf, 30, 35, 40, 45, Inf),
      labels = c("<30", "30-35", "35-40", "40-45", "≥45"),
      right = FALSE # Includes lower bound: [30, 35)
    ),
    # Drop matrix attributes created by scale()
    lat_centred = as.numeric(scale(lat, scale = FALSE)),
    lon_centred = as.numeric(scale(lon, scale = FALSE))
  )

## Predicted Daymet
predicted_daymet <- bind_rows(predicted_daymet1, predicted_daymet2) %>%
  distinct(RouteNumber, year, .keep_all = TRUE)

obs_pred_dfc_daymet <- observed %>%
  left_join(predicted_daymet, join_by(
    RouteNumber == RouteNumber,
    SurveyYear == year
  )) %>%
  select(-intensity) %>%
  janitor::clean_names() %>%
  rename(
    observed_dfc = "onset_doy",
    predicted_dfc = "first"
  ) %>%
  mutate(
    # Create ordered latitude bins safely
    lat_bin = cut(
      lat,
      breaks = c(-Inf, 30, 35, 40, 45, Inf),
      labels = c("<30", "30-35", "35-40", "40-45", "≥45"),
      right = FALSE # Includes lower bound: [30, 35)
    ),
    # Drop matrix attributes created by scale()
    lat_centred = as.numeric(scale(lat, scale = FALSE)),
    lon_centred = as.numeric(scale(lon, scale = FALSE))
  )

# Difference in distribution of predicted DFC in ERA5 Land vs. Daymet
hist(obs_pred_dfc$predicted_dfc)
hist(obs_pred_dfc_daymet$predicted_dfc)
## Looks similar but Daymet predicted DFCs are fewer in the middle and are skewed slightly earlier

sd(obs_pred_dfc$predicted_dfc)
sd(obs_pred_dfc_daymet$predicted_dfc)
## spread of data is very similar - almost identical ERA5 - 22.686; Daymet 22.330

# Build linear models ----------------------------------------------------------
# Excluded 'route' random effect; high site uniqueness (74/90, ~82%) provides insufficient repeated measures to estimate group-level variance

# Model with just TS3 as predictor for observed calling
mod <- lm(observed_dfc ~ predicted_dfc, data = obs_pred_dfc)
summary(mod)
model_offset <- lm(observed_dfc ~ predicted_dfc + offset(1 * predicted_dfc), data = obs_pred_dfc)
summary(model_offset)
vis_mod <- ggeffects::ggpredict(mod, terms = c("predicted_dfc"))

# Plot
fig3 <- ggplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "#FFFFFF") +
  # Raw data points
  geom_point(
    data = obs_pred_dfc,
    mapping = aes(x = predicted_dfc, y = observed_dfc, color = lat_bin, shape = lat_bin),
    size = 2,
    alpha = 0.7
  ) +
  # Model 95% Confidence Interval
  geom_ribbon(
    data = vis_mod,
    mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.3,
    fill = "#FFFFFF"
  ) +
  # Model Fitted Line
  geom_line(
    data = vis_mod,
    mapping = aes(x = x, y = predicted),
    colour = "#FFFFFF"
  ) +
  theme_classic() +
  theme(
    text = element_text(colour = "#FFFFFF"),
    axis.text = element_text(colour = "#FFFFFF"),
    axis.line = element_line(colour = "#FFFFFF"),
    axis.ticks = element_line(colour = "#FFFFFF"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  labs(x = "Predicted DFC Based on TS3", y = "Observed DFC", colour = "Latitude", shape = "Latitude") +
  scale_colour_viridis_d(begin = 0.6, end = 0.9, direction = -1)

ggsave("figs/figure3.png", plot = fig3, width = 6.5, height = 4, dpi = 400, bg = "transparent")

fig3a <- ggplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  # Raw data points
  geom_point(
    data = obs_pred_dfc,
    mapping = aes(x = predicted_dfc, y = observed_dfc, color = lat_bin, shape = lat_bin),
    size = 2,
    alpha = 0.7
  ) +
  # Model 95% Confidence Interval
  geom_ribbon(
    data = vis_mod,
    mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.3
  ) +
  # Model Fitted Line
  geom_line(
    data = vis_mod,
    mapping = aes(x = x, y = predicted)
  ) +
  theme_classic() +
  labs(x = "Predicted DFC Based on TS3", y = "Observed DFC", colour = "Latitude", shape = "Latitude") +
  scale_colour_viridis_d(begin = 0.6, end = 0.9, direction = -1)

ggsave("figs/figure3a.png", plot = fig3a, width = 6.5, height = 4, dpi = 400, bg = "transparent")

# Model with survey year as predictor
mod2 <- lm(observed_dfc ~ predicted_dfc + survey_year, data = obs_pred_dfc)
summary(mod2)
vis_mod2 <- ggeffects::ggpredict(mod2, terms = c("survey_year"))

ggplot() +
  # Raw data points over time
  geom_point(
    data = obs_pred_dfc,
    mapping = aes(x = survey_year, y = observed_dfc),
    alpha = 0.2
  ) +
  # Model 95% CI band for year effect
  geom_ribbon(
    data = vis_mod2,
    mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.3
  ) +
  # Model line across years (flat line = no systematic year effect)
  geom_line(
    data = vis_mod2,
    mapping = aes(x = x, y = predicted),
    linewidth = 1
  ) +
  theme_bw() +
  labs(
    x = "Survey Year",
    y = "Observed DFC",
  )

# Model with survey_year as a random effect
mod_year_random <- lmer(observed_dfc ~ predicted_dfc + (1 | survey_year), data = obs_pred_dfc)
summary(mod_year_random)
performance::r2(mod_year_random)
## estimate is stronger when year is included as a random effect
## conditional R2 0.321, marginal R2 0.219

# Plot
year_re <- broom.mixed::tidy(mod_year_random,
  effects   = "ran_vals",
  conf.int  = TRUE
) %>%
  left_join(
    obs_pred_dfc %>%
      count(survey_year) %>%
      mutate(level = as.character(survey_year)),
    by = "level"
  ) %>%
  mutate(
    year_label = paste0(level, "  (n = ", n, ")"),
    year_label = factor(year_label, levels = rev(sort(unique(year_label))))
  )

ggplot(year_re, aes(x = estimate, y = year_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
    size = 0.4
  ) +
  labs(
    x = "Deviation from mean calling onset (days)",
    y = NULL,
    title = "Interannual variation in calling onset",
    subtitle = "Year-level random intercepts, sample size shown per year"
  ) +
  theme_bw()

# model with lat, lon, and survey_year as a random effect
mod_year_random2 <- lmer(observed_dfc ~ predicted_dfc + lat_centred + lon_centred + (1 | survey_year), data = obs_pred_dfc)
summary(mod_year_random2)

## when you include lat and lon with survey year as a random effect, you get a singular fit, meaning that survey_year as a fixed effect does not include any information beyond what is currently included in the model, or there are not enough values per year to separate variance

# model with lat and lon as predictors
mod3 <- lm(observed_dfc ~ predicted_dfc + survey_year + lat_centred + lon_centred, data = obs_pred_dfc)
summary(mod3)
# note that the year effect from mod2 has disappeared - which means this confirms a sampling composition effect rather than an actual systematic year change

# drop survey_year
mod4 <- lm(observed_dfc ~ predicted_dfc + lat_centred + lon_centred, data = obs_pred_dfc)
summary(mod4)

# Plot
vis_mod4 <- ggeffects::ggpredict(mod4, terms = c("predicted_dfc"))

fig4 <- ggplot() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  # Raw data points over time
  geom_point(
    data = obs_pred_dfc,
    mapping = aes(x = predicted_dfc, y = observed_dfc),
    alpha = 0.2
  ) +
  # Model 95% CI band for year effect
  geom_ribbon(
    data = vis_mod4,
    mapping = aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.3
  ) +
  # Model line across years (flat line = no systematic year effect)
  geom_line(
    data = vis_mod4,
    mapping = aes(x = x, y = predicted),
    linewidth = 1
  ) +
  theme_bw() +
  labs(
    x = "Predicted DFC Based on TS3 - with Constant Lat and Lon",
    y = "Observed DFC",
  )

ggsave("figs/figure4.png", plot = fig4, width = 6.5, height = 4, dpi = 400)

# model with only lat and lon as predictors
mod5 <- lm(observed_dfc ~ lat_centred + lon_centred, data = obs_pred_dfc)
summary(mod5)

# Check model parsimony with AIC
AIC(mod, mod2, mod3, mod4, mod5)

# AIC confirms that geography + thermal sum is best model - check to see how much variance is accounted for with ANOVA
anova(mod5, mod4)
## significant amount of variance accounted for with the inclusion of predicted_dfc - makes sense since slope was significant in model output

# Mod4 is more parsimonious, but the Estimate for predicted_dfc actually shows a negative effect instead of a positive one. This may be due to (1) the lower latitude sites that were shown to have a floor effect - DFC is reached very early after February 1st (start of the thermal sum calculation window), (2) influential observations may be skewing results, (3) there is high collinearity between observations, (4), the ERA5 Land hourly resolution is too coarse

# (1) Do a sensitivity analysis to see if removing the sites at lat_bin 30 changes the amount of variance explained by the thermal sum (TS3)

mid_obs_pred_dfc <- obs_pred_dfc %>%
  filter(lat_bin != "30-35") %>%
  mutate(
    lat_centred = scale(lat, scale = FALSE),
    lon_centred = scale(lon, scale = FALSE)
  )

mod4_sensitivity <- lm(observed_dfc ~ predicted_dfc + lat_centred + lon_centred, data = mid_obs_pred_dfc)
summary(mod4_sensitivity)
car::Anova(mod4, mod4_sensitivity)
## Removing the lat-bin 30 sites did not move the estimate much. Rejected.

# (2) Look for influential observations with Cooks distance
plot(mod4, which = 5)
max(cooks.distance(mod4))
mean(cooks.distance(mod4))
sum(cooks.distance(mod4) > 4 / nrow(obs_pred_dfc))
## 5 observations are flagged by sum(cooks D) > 4 / n
## max leverage ~0.085 vs mean ~0.009, which is below 1

# Check stricter threshold
cooks_vals <- cooks.distance(mod4)
plot(cooks_vals,
  type = "h", main = "Cook's Distance",
  ylab = "Cook's Distance", xlab = "Observation Index", lwd = 2
)
n <- nrow(obs_pred_dfc)
threshold <- 4 / n
abline(h = threshold, lty = 2, lwd = 2)
## Observations 3, 15, 89, 90, 91 are flagged - see if they change the model
# Rerun the model excluding rows
mod4_clean <- lm(observed_dfc ~ predicted_dfc + lat_centred + lon_centred, data = obs_pred_dfc[-c(3, 15, 89, 90, 91), ])

# Compare the old and new results
summary(mod4)
summary(mod4_clean)
## More negative, influential points not driving TS3 negative coefficient

# (3) Assess for collinearity between predictors
car::vif(mod4)
summary(lm(predicted_dfc ~ lat_centred + lon_centred, data = obs_pred_dfc))
## Interesting that the variance inflation factor doesn't seem to be too high, but still around 75% of variance in predicted_dfc can be reconstructed by lat and lon.

# (4) Climate product resolution effects
# Use predicted DFC calculated from Daymet interpolated climate dataset with ~1km resolution to see whether the TS3 index estimate changes
# how much of each predictor is geography?
summary(lm(predicted_dfc ~ lat_centred + lon_centred, data = obs_pred_dfc_daymet))$r.squared
## Similar to the ERA5 Land predicted dfc, ~73% of the predictor is explained by geography

# does the finer resolution change the effectiveness of the thermal index
summary(lm(observed_dfc ~ predicted_dfc + lat_centred + lon_centred, data = obs_pred_dfc_daymet))
## predicted_dfc still produces a negative estimate, finer resolution likely not the cause if routes are farther apart, ERA5 still produces the same result. Rejected

# Model geography + TS3 with an interaction term
mod6 <- lm(observed_dfc ~ predicted_dfc * lat_centred + lon_centred, data = obs_pred_dfc)
summary(mod6)
anova(mod4, mod6)
## looks like interaction is not significant, meaning that TS3's effect does not seem to differ when changing latitudes - this may just be an artifact of sampling - worth discussing!

# Generate summary tables -------------------------------------------------------
# Summarize selection criteria for best-fit among models
models <- list(
  "TS3" = mod,
  "TS3 + year" = mod2,
  "Latitude + longitude" = mod5,
  "TS3 + latitude + longitude" = mod4,
  "TS3 × latitude + longitude" = mod6
)

tibble(
  Model = names(models),
  k     = map_dbl(models, ~ length(coef(.x))),
  R2    = map_dbl(models, ~ summary(.x)$r.squared),
  RSE   = map_dbl(models, ~ summary(.x)$sigma),
  AIC   = map_dbl(models, AIC)
) %>%
  mutate(dAIC = AIC - min(AIC)) %>%
  select(-AIC) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

# Summarize model 4 (best-fit)
broom::tidy(mod4) %>%
  mutate(
    term = recode(term,
      "(Intercept)" = "Intercept",
      predicted_dfc = "Predicted DFC (TS3)",
      lat_centred = "Latitude (days per °N)",
      lon_centred = "Longitude (days per °E)"
    ),
    across(c(estimate, std.error), ~ round(.x, 3)),
    p.value = ifelse(p.value < 0.001, "< 0.001", sprintf("%.3f", p.value))
  ) %>%
  select(Term = term, Estimate = estimate, SE = std.error, p = p.value)