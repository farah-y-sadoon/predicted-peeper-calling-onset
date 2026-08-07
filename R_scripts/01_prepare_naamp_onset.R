# Prepare North American Amphibian Monitoring Program (NAAMP) Data for Model Validation
# Get nearest NAAMP site to Lovett's study site (2013)
# Use nearest NAAMP site to download predicted DFC values from GEE

# Load Packages ----
library(tidyverse)
library(dplyr)
library(lubridate)
library(sf)

# Data import -----------------------------------------------------------------

runs <- read_csv("data_raw/Runs.csv")
stops <- read_csv("data_raw/Stops.csv")
counts <- read_csv("data_raw/Counts.csv")
species <- read_csv("data_raw/Species.csv")
coords <- read_csv("data_raw/Coordinates.csv")

# Get onset of calling --------------------------------------------------------
# Add julian date
runs2 <- runs %>%
  mutate(
    SurveyDate = mdy(SurveyDate),
    DOY = yday(SurveyDate)
  )

# join to each run
counts2 <- counts %>%
  left_join(
    runs2 %>%
      select(RunID, RouteNumber, SurveyYear, DOY, RunNumber, QuizScore, UnifiedProtocol),
    by = "RunID"
  )

# Test ------------------------------------------------------------------------

# Calculate onset for every species, route, year, and threshold
onset_route <- counts2 %>%
  group_by(Species, RouteNumber, SurveyYear) %>%
  summarise(
    onset_1 = if(any(CallingIndex >= 1))
      min(DOY[CallingIndex >= 1]) else NA_real_,
    onset_2 = if(any(CallingIndex >= 2))
      min(DOY[CallingIndex >= 2]) else NA_real_,
    onset_3 = if(any(CallingIndex >= 3))
      min(DOY[CallingIndex >= 3]) else NA_real_,
    .groups = "drop"
  )
thresholds <- c(1, 2, 3)

onset_route <- purrr::map_dfr(thresholds, function(thresh){
  
  counts2 %>%
    filter(CallingIndex >= thresh) %>%
    group_by(Species, RouteNumber, SurveyYear) %>%
    summarise(
      onset_DOY = min(DOY),
      .groups = "drop"
    ) %>%
    mutate(threshold = thresh)
  
})

head(onset_route)

# Quality control -------------------------------------------------------------
# if species werent detected, they were not written down. so if the run happened and something called but not the other species, then it shoudl be 0

# all surveys that occurred
survey_grid <- runs2 %>%
  distinct(RouteNumber, SurveyYear, RunID, DOY) %>%
  crossing(
    Species = unique(counts2$Species)
  )
survey_grid <- runs2 %>%
  select(
    RunID,
    RouteNumber,
    SurveyYear,
    DOY,
    RunNumber,
    QuizScore,
    UnifiedProtocol
  ) %>%
  distinct() %>%
  crossing(
    Species = unique(counts2$Species)
  )

# observed calling intensity
calling <- counts2 %>%
  group_by(Species, RouteNumber, SurveyYear, RunID, DOY) %>%
  summarise(
    MaxCalling = max(CallingIndex, na.rm = TRUE),
    .groups = "drop"
  )

# add non-detections as zero  
# summarize by route run (collapse stops)

route_run <- survey_grid %>%
  left_join(
    calling,
    by = c("Species", "RouteNumber", "SurveyYear", "RunID", "DOY")
  ) %>%
  mutate(
    MaxCalling = replace_na(MaxCalling, 0)
  )

route_run %>% count(MaxCalling) #working because now we see a bunch of 0s

# The onset is the first survey where calling reaches a threshold
# provided that at least one previous survey had a lower calling index
# and has not already exceeded it
# For onset of intensity 1: first transition < 1 to >/= 1
# For onset of intensity 2: first transition < 2 to >/=2
# For onset of intensity 3: first transition < 3 >= 3
## this was why we needed to put the 0s in

# and track the number of days between runs for later QC
onset_route <- route_run %>%
  arrange(Species, RouteNumber, SurveyYear, DOY) %>%
  group_by(Species, RouteNumber, SurveyYear) %>%
  mutate(
    prev_call = lag(MaxCalling),
    days_since_prev = DOY - lag(DOY)
  ) %>%
  summarise(
    
    onset_1 = {
      x <- which(MaxCalling >= 1 & prev_call < 1)
      if(length(x)>0) DOY[min(x)] else NA_real_
    },
    
    gap_1 = {
      x <- which(MaxCalling >= 1 & prev_call < 1)
      if(length(x)>0) days_since_prev[min(x)] else NA_real_
    },
    
    quiz_1 = {
      x <- which(MaxCalling >= 1 & prev_call < 1)
      if(length(x)>0) QuizScore[min(x)] else NA
    },
    
    protocol_1 = {
      x <- which(MaxCalling >= 1 & prev_call < 1)
      if(length(x)>0) UnifiedProtocol[min(x)] else NA
    },
    
    
    onset_2 = {
      x <- which(MaxCalling >= 2 & prev_call < 2)
      if(length(x)>0) DOY[min(x)] else NA_real_
    },
    
    gap_2 = {
      x <- which(MaxCalling >= 2 & prev_call < 2)
      if(length(x)>0) days_since_prev[min(x)] else NA_real_
    },
    
    quiz_2 = {
      x <- which(MaxCalling >= 2 & prev_call < 2)
      if(length(x)>0) QuizScore[min(x)] else NA
    },
    
    protocol_2 = {
      x <- which(MaxCalling >= 2 & prev_call < 2)
      if(length(x)>0) UnifiedProtocol[min(x)] else NA
    },
    
    
    onset_3 = {
      x <- which(MaxCalling >= 3 & prev_call < 3)
      if(length(x)>0) DOY[min(x)] else NA_real_
    },
    
    gap_3 = {
      x <- which(MaxCalling >= 3 & prev_call < 3)
      if(length(x)>0) days_since_prev[min(x)] else NA_real_
    },
    
    quiz_3 = {
      x <- which(MaxCalling >= 3 & prev_call < 3)
      if(length(x)>0) QuizScore[min(x)] else NA
    },
    
    protocol_3 = {
      x <- which(MaxCalling >= 3 & prev_call < 3)
      if(length(x)>0) UnifiedProtocol[min(x)] else NA
    },
    
    .groups="drop"
  )
head(onset_route)

onset_route %>%
  summarise(
    n1 = sum(!is.na(onset_1)),
    n2 = sum(!is.na(onset_2)),
    n3 = sum(!is.na(onset_3))
  )

# identify where gaps between surveys were larger than 30 days, do not include
# in final tally of onset

onset_route_qc <- onset_route %>%
  mutate(
    onset_1_raw = onset_1,
    onset_2_raw = onset_2,
    onset_3_raw = onset_3,
    
    onset_1 = ifelse(gap_1 <= 30, onset_1, NA_real_),
    onset_2 = ifelse(gap_2 <= 30, onset_2, NA_real_),
    onset_3 = ifelse(gap_3 <= 30, onset_3, NA_real_)
  )

onset_route %>%
  summarise(
    onset1 = sum(!is.na(onset_1)),
    onset2 = sum(!is.na(onset_2)),
    onset3 = sum(!is.na(onset_3))
  )

onset_route_qc %>%
  summarise(
    onset1 = sum(!is.na(onset_1)),
    onset2 = sum(!is.na(onset_2)),
    onset3 = sum(!is.na(onset_3))
  )

# Organize for downstream analysis --------------------------------------------
# tidy for plotting and further analysis
onset_long <- onset_route_qc %>%
  pivot_longer(
    cols = starts_with("onset_"),
    names_to = "intensity",
    values_to = "onset_DOY"
  )

# select only peepers
peepers <- onset_long %>%
  filter(Species == "Pseudacris crucifer") %>% 
  select(RouteNumber, SurveyYear, intensity, onset_DOY) %>% 
  filter(intensity %in% c('onset_1', 'onset_2', 'onset_3'), 
         !is.na(onset_DOY))

# get centroid for each route's coordinates 
route_coords <- coords %>% 
  group_by(RouteNumber) %>% 
  summarise(lat = mean(lat), 
            lon = mean(lon),
            .groups = "drop")

# join peepers with route coordinates
peepers <- peepers %>% 
  left_join(route_coords, by = c("RouteNumber"))

# filter out NA values
peepers <- peepers %>% 
  filter(!is.na(lat), 
         !is.na(lon))

# define Lovett 2013 coordinates for study site and create spatial feature
lovett_site <- list(lon  = -(73 + (44.90 / 60)),
                    lat  = 41 + (51.01 / 60))

target_geom <- st_sfc(st_point(c(lovett_site$lon, lovett_site$lat)), crs = 4326)

# extract only coordinates for peeper routes, create spatial features and find the point closest to Lovett's site
peeper_route_coords <- peepers %>% 
  select(RouteNumber, lat, lon) %>% 
  distinct()

peeper_routes_sf <- st_as_sf(peeper_route_coords, coords = c("lon", "lat"), crs = 4326)

# find NAAMP route closest to Lovett's site
closest_index <- st_nearest_feature(st_sfc(target_geom), peeper_routes_sf)
closest_route <- peeper_routes_sf[closest_index, ] # Route 610323 POINT (-73.77666 41.81619)

# calculate distance between Lovett's site and closest point
closest_route_dist <- st_distance(target_geom, closest_route) # 4447.991m (4.448km)