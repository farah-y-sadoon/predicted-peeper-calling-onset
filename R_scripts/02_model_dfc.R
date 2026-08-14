# ==============================================================================
# Modelling NAAMP Peeper Calling Onset with Lovett (2013) Predicted DFC
# ==============================================================================

# Load Packages ----------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(lubridate)
library(lme4)

# Data import ------------------------------------------------------------------
## NAAMP Observations
observed <- read_csv("data_processed/naamp_peepers_onset.csv")

## Predicted DFC from ERA5 Land (~9km resolution, ~81km2 )
predicted1 <- read_csv("data_processed/naamp_route_predicted_dfc_1997_2005.csv")
predicted2 <- read_csv("data_processed/naamp_route_predicted_dfc_2006_2015.csv")

## Predicted DFC from Daymet (~1km resolution, 1km2)
predicted_daymet1 <- read_csv("data_processed/naamp_route_predicted_dfc_daymet_1997_2005.csv")
predicted_daymet2 <- read_csv("data_processed/naamp_route_predicted_dfc_daymet_2006_2015.csv")

# Data prep --------------------------------------------------------------------
## Predicted ERA5 Land
predicted <- bind_rows(predicted1, predicted2) %>%
  distinct(RouteNumber, year, .keep_all = TRUE)

obs_pred_dfc <- observed %>% 
  left_join(predicted, join_by(RouteNumber == RouteNumber, 
                               SurveyYear == year)) %>% 
  select(-intensity) %>% 
  janitor::clean_names() %>% 
  rename(observed_dfc = "onset_doy", 
         predicted_dfc = "first") %>% 
  mutate(lat_bin = case_when(lat >= 45 ~ "45", 
                             lat >= 40 ~ "40", 
                             lat >= 35 ~ "35", 
                             lat >= 30 ~ "30"), 
         lat_centred = scale(lat, scale = FALSE), 
         lon_centred = scale(lon, scale = FALSE))

## Predicted Daymet
predicted_daymet <- bind_rows(predicted_daymet1, predicted_daymet2) %>% 
  distinct(RouteNumber, year, .keep_all = TRUE)

obs_pred_dfc_daymet <- observed %>% 
  left_join(predicted_daymet, join_by(RouteNumber == RouteNumber, 
                                       SurveyYear == year)) %>% 
  select(-intensity) %>% 
  janitor::clean_names() %>% 
  rename(observed_dfc = "onset_doy", 
         predicted_dfc = "first") %>% 
  mutate(lat_bin = case_when(lat >= 45 ~ "45", 
                             lat >= 40 ~ "40", 
                             lat >= 35 ~ "35", 
                             lat >= 30 ~ "30"), 
         lat_centred = scale(lat, scale = FALSE), 
         lon_centred = scale(lon, scale = FALSE))

# Difference in distribution of predicted DFC in ERA5 Land vs. Daymet
hist(obs_pred_dfc$predicted_dfc)
hist(obs_pred_dfc_daymet$predicted_dfc)

sd(obs_pred_dfc$predicted_dfc)
sd(obs_pred_dfc_daymet$predicted_dfc)
## Looks like 

# Build linear model -----------------------------------------------------------
# route was not included as a random effect because 74/90 sites in the model are unique (~82%)
# meaning that pseudoreplication is not a major issue in the data

# model with just predicted as predictor for observed calling
mod <- lm(observed_dfc ~ predicted_dfc, data = obs_pred_dfc)
summary(mod)

ggplot(obs_pred_dfc, aes(predicted_dfc, observed_dfc)) +
  geom_point(size = 3, aes(predicted_dfc, observed_dfc, colour = lat_bin)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm") +
  theme_classic()
### UPDATE: calculate confidence intervals and plot them manually like in Megan's script

# model with survey year as fixed effect
mod2 <- lm(observed_dfc ~ predicted_dfc + survey_year, data = obs_pred_dfc)
summary(mod2)

# model with lat and lon as fixed effects
mod3 <- lm(observed_dfc ~ predicted_dfc + survey_year + lat_centred + lon_centred, data = obs_pred_dfc)
summary(mod3)

ggplot(obs_pred_dfc, aes(lat, observed_dfc)) + 
  geom_point(size = 2, alpha = 0.5) + 
  geom_smooth(method = "lm") + 
  theme_classic()
### UPDATE: calculate confidence intervals and plot them manually like in Megan's script
# note that the year effect from mod2 has disappeared - which means this is likely a sampling composition effect rather than an actual year change

# drop survey_year 
mod4 <- lm(observed_dfc ~ predicted_dfc + lat_centred + lon_centred, data = obs_pred_dfc)
summary(mod4)

# model with only lat and lon as fixed effects
mod5 <- lm(observed_dfc ~ lat_centred + lon_centred, data = obs_pred_dfc)
summary(mod5)

# Check model parsimony with AIC
AIC(mod, mod2, mod3, mod4, mod5)

# AIC confirms that geography + thermal sum is best model - check to see how much variance is accounted for with ANOVA

anova(mod5, mod4)

# Mod4 is more parsimonious, but the Estimate for predicted_dfc actually shows a negative effect instead of a positive one. This may be due to (1) the lower latitude sites that were shown to have a floor effect - DFC is reached very early after February 1st (start of the thermal sum calculation window), (2), (3), (4)

# (1) Do a sensitivity analysis to see if removing the sites at lat_bin 30 changes the amount of variance explained by the thermal sum (TS3)

mid_obs_pred_dfc <- obs_pred_dfc %>% 
  filter(lat_bin != 30) %>% 
  mutate(lat_centred = scale(lat, scale = FALSE), 
         lon_centred = scale(lon, scale = FALSE))

summary(lm(observed_dfc ~ predicted_dfc + lat_centred + lon_centred, data = mid_obs_pred_dfc))

## Removing the lat-bin 30 sites did not move the estimate much - could be due to highly collinear relationship between the thermal sum value and latitude - so TS3 does predict, but does not add much value beyond what geography already tells us

# (2) Look for influential observations with Cooks distance
plot(mod4, which = 5)
max(cooks.distance(mod4))
mean(cooks.distance(mod4))
sum(cooks.distance(mod4) > 4 / nrow(obs_pred_dfc))
# No points near Cook's D thresholds; max leverage ~0.085 vs mean ~0.009. Rejected.

# (3) Assess for collinearity between predictors 
car::vif(mod4)
summary(lm(predicted_dfc ~ lat_centred + lon_centred, data = obs_pred_dfc))
## Interesting that the variance inflation factor doesn't seem to be too high, but still around 75% of variance in predicted_dfc can be reconstructed by lat and lon

# (4) Climate product resolution effects
# Use predicted DFC calculated from Daymet interpolated climate dataset with ~1km resolution to see whether the TS3 index estimate changes
# how much of each predictor is geography?
summary(lm(predicted_dfc ~ lat_centred + lon_centred, data = obs_pred_dfc_daymet))$r.squared
## Similar to the ERA5 Land predicted dfc, ~73% of the predictor is explained by geography 

# does the finer resolution change the effectiveness of the thermal index
summary(lm(observed_dfc ~ predicted_dfc + lat_centred + lon_centred, data = obs_pred_dfc_daymet))
## predicted_dfc still produces a negative estimate, finer resolution likely not the cause if routes are farther apart, ERA5 still produces the same result
