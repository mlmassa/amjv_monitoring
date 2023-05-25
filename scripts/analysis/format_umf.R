# Format data for occupancy analysis

# Setup -------------------------------------------------------------------

library(tidyverse)
library(suncalc)
library(unmarked)

theme_set(
  theme_minimal(base_size = 11) +
    theme(strip.background = element_blank()))

# Set parameters
focal_spp <- 
  c("WOTH", "WEWA", "SCTA", "YBCU", # focus for abund models
    "CERW", "CAWA", "BBCU", "RUGR", "KEWA") # Not abund modelable

# Import data -------------------------------------------------------------

# NOTE: 5 visits occurred after instructed cutoff time, all 2021 Fed points

counts <- 
  read_rds("data/raw/amjv_data.rds")$count_data

visits <- 
  read_rds("data/raw/amjv_data.rds")$visits |> 
  mutate(year = year(date))

points <- 
  read_rds("data/raw/amjv_data.rds")$points

plots <- 
  read_rds("data/raw/amjv_data.rds")$plots

veg <- 
  read_rds("data/raw/amjv_data.rds")$veg_data

# Subset to relevant variables
timefix <-
  visits |> 
  select(point_id, visit_id, date, start_time) |> 
  # Add coordinates of point
  left_join(points |> select(point_id, lat, long), by = "point_id") |> 
  # Add sunrise time
  left_join(
    getSunlightTimes(
      data = 
        visits |> 
        select(point_id, date, start_time) |> 
        mutate(date = as.Date(date)) |> 
        left_join(
          points |> dplyr::select(point_id, lat, long), 
          by = "point_id") |> 
          select(date, lat, lon = long),
      keep = "sunrise",
      tz = "America/New_York") |> 
      as_tibble(),
    by = c("date", "lat", "long" = "lon"))

date(timefix$start_time) <- timefix$date
tz(timefix$start_time) <- "America/New_York"

visits <-
  visits |> 
  left_join(
    timefix |> 
      # Calculate time since sunrise
      mutate(
        diff = difftime(start_time, sunrise, units = "mins"),
        # Convert to numeric, "minutes after sunrise"
        start_sun = as.numeric(diff, units = "mins")) %>% 
      # Select only calculation, then append to visits
      select(visit_id, start_sun),
    by = "visit_id")

# Generate untidy dataset:
birds <-
  counts |> 
  left_join(visits, by = "visit_id") |> 
  left_join(points, by = "point_id") |> 
  left_join(plots, by = "plot_name")

# Generate detection table ------------------------------------------------

# Create annual point visit table
annual_points <-
  visits |> 
  # Get unique points for each year
  distinct(point_id, year) |> 
  # Mark these as visited
  mutate(visited = 1) |> 
  # Get full table of points per year
  complete(year, nesting(point_id), fill = list(visited = 0)) |> 
  arrange(point_id, year)

# Generate detection history
detections <-
  birds |> 
  # No flyover detections: birds in habitat only
  filter(flyover == 0) |> 
  # For each species at each pointxyearxvisit...
  group_by(point_id, year, visit, species) |> 
  # Get abundance
  summarize(abund = n(), .groups = "drop") |> 
  # Add empty combinations not present
  right_join(
    birds |> expand(point_id, year, visit, species),
    by = c("point_id", "year", "visit", "species")) |> 
  # Add whether pointxyear combo was visited
  left_join(annual_points) |>
  # Add detection (1/0 if visited, keep NA if not)
  mutate(
    abund =
      if_else(
        is.na(abund) & visited == 1, 0,
        abund),
    det = 
      case_when(
        is.na(abund) & visited == 1 ~ 0,
        abund == 0 ~ 0,
        abund > 0 ~ 1)) |> 
  # Retain only focal species
  filter(species %in% focal_spp)

# Plot abundances (zero-inflated)
detections |> 
  filter(!is.na(abund)) |> 
  ggplot() + 
  geom_bar(aes(x = factor(abund))) + 
  facet_wrap(vars(species)) +
  labs(x = "Abundance (birds/survey)")

ggsave(
  "output/plots/focal_abund.png",
  width = 6, height = 5, units = "in", dpi = 300)

# Plot occupancy
detections |> 
  filter(!is.na(abund)) |> 
  ggplot() + 
  geom_bar(aes(x = factor(det))) + 
  facet_wrap(vars(species)) +
  labs(x = "Detected")

ggsave(
  "output/plots/focal_det.png",
  width = 6, height = 5, units = "in", dpi = 300)

occupancy <-
  detections |> 
  arrange(point_id, year, visit)  |> 
  # Cut irrelevant variables
  select(point_id, year, visit, species, det) |> 
  # Make wide data with column for each year-visit detection
  # detection is "y"
  pivot_wider(
    names_prefix = 'y_', 
    names_from = visit,
    values_from = det,
    names_sep = '_') |> 
  unite("point_year", c(point_id, year), remove = T) |> 
  arrange(point_year, species)

abundance <-
  detections |> 
  arrange(point_id, year, visit)  |> 
  # Cut irrelevant variables
  select(point_id, year, visit, species, abund) |> 
  # Make wide data with column for each year-visit detection
  # detection is "y"
  pivot_wider(
    names_prefix = 'y_', 
    names_from = visit,
    values_from = abund,
    names_sep = '_') |> 
  unite("point_year", c(point_id, year), remove = T) |> 
  arrange(point_year, species)

# Generate observation covariates -----------------------------------------

# Format:
# siteCovs: single (n_sitexyear) row x n_covariates column DF 
# obsCovs: list of (n_sitexyear) row x n_visits column DFs

## Observer (factor) ----

ggplot(data = visits) +
  geom_bar(aes(x = observer)) +
  facet_wrap(vars(year), nrow = 2) +
  labs(
    y = "Surveys",
    x = "Observer") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/observer.png",
  width = 3.5, height = 3.5, units = "in", dpi = 300)

observer <-
  visits |> 
  select(point_id, visit, year, observer) |> 
  right_join(
    birds |> expand(point_id, year, visit),
    by = c("point_id", "year", "visit")) |> 
  mutate(observer = as.factor(observer)) |> 
  pivot_wider(
    names_prefix = 'observer_',
    names_from = c(visit),
    values_from = observer,
    names_sep = '_') |> 
  unite("point_year", c(point_id, year), remove = T) |> 
  arrange(point_year)

## Wind (factor) ----

ggplot(data = visits) +
  geom_bar(aes(x = factor(wind))) +
  facet_wrap(vars(year), nrow = 2) +
  labs(
    y = "Surveys",
    x = "Wind (Beaufort scale)") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/wind.png",
  width = 3.5, height = 3.5, units = "in", dpi = 300)

wind <-
  visits |> 
  select(point_id, visit, year, wind) |> 
  right_join(
    birds |> expand(point_id, year, visit),
    by = c("point_id", "year", "visit")) |> 
  mutate(wind = as.factor(wind)) |> 
  pivot_wider(
    names_prefix = 'wind_',
    names_from = c(visit),
    values_from = wind,
    names_sep = '_') |> 
  unite("point_year", c(point_id, year), remove = T) |> 
  arrange(point_year)

## Temperature (c) (2022 only) ----

ggplot(data = visits |> filter(!is.na(temp))) +
  geom_histogram(aes(x = temp)) +
  facet_wrap(vars(year), nrow = 2) +
  labs(
    y = "Surveys",
    x = "Temperature (C)") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/temp.png",
  width = 3.5, height = 3.5, units = "in", dpi = 300)

temp <-
  visits |> 
  select(point_id, visit, year, temp) |> 
  right_join(
    birds |> expand(point_id, year, visit),
    by = c("point_id", "year", "visit")) |> 
  pivot_wider(
    names_prefix = 'temp_',
    names_from = c(visit),
    values_from = temp,
    names_sep = '_') |> 
  unite("point_year", c(point_id, year), remove = T) |> 
  arrange(point_year)

## Precip (factor) ----

ggplot(data = visits) +
  geom_bar(aes(x = precip)) +
  facet_wrap(vars(year), nrow = 2) +
  labs(
    y = "Surveys",
    x = "Precipitation") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/precip.png",
  width = 3.5, height = 3.5, units = "in", dpi = 300)

precip <-
  visits |> 
  select(point_id, visit, year, precip) |> 
  mutate(precip = factor(precip)) |> 
  right_join(
    birds |> expand(point_id, year, visit),
    by = c("point_id", "year", "visit")) |> 
  pivot_wider(
    names_prefix = 'precip_',
    names_from = c(visit),
    values_from = precip,
    names_sep = '_') |> 
  unite("point_year", c(point_id, year), remove = T) |> 
  arrange(point_year)

## Cloud ----

ggplot(data = visits) +
  geom_histogram(aes(x = cloud), binwidth = 10) +
  facet_wrap(vars(year), nrow = 2) +
  labs(
    y = "Surveys",
    x = "Cloud cover (%)") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/cloud.png",
  width = 3.5, height = 3.5, units = "in", dpi = 300)

cloud <-
  visits |> 
  select(point_id, visit, year, cloud) |> 
  right_join(
    birds |> expand(point_id, year, visit),
    by = c("point_id", "year", "visit")) |> 
  pivot_wider(
    names_prefix = 'cloud_',
    names_from = c(visit),
    values_from = cloud,
    names_sep = '_') |> 
  unite("point_year", c(point_id, year), remove = T) |> 
  arrange(point_year)

## Minutes since sunrise ----

ggplot(data = visits) +
  geom_histogram(aes(x = start_sun), binwidth = 10) +
  # geom_vline(xintercept = (4*60)) +
  facet_wrap(vars(year), nrow = 2) +
  scale_x_continuous(breaks = c(0, 60, 120, 180, 240, 300)) +
  labs(
    y = "Surveys",
    x = "Start time (minutes after sunrise)") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/start_sun.png",
  width = 3.5, height = 3.5, units = "in", dpi = 300)

start_sun <-
  visits |> 
  select(point_id, visit, year, start_sun) |> 
  right_join(
    birds |> expand(point_id, year, visit),
    by = c("point_id", "year", "visit")) |> 
  pivot_wider(
    names_prefix = 'start_sun_',
    names_from = c(visit),
    values_from = start_sun,
    names_sep = '_') |> 
  unite("point_year", c(point_id, year), remove = T) |> 
  arrange(point_year)

## Day of year ----
ggplot(
  data = visits |> 
    mutate(doy = as.Date(yday(date) - 1, origin = "2023-01-01"))) +
  scale_x_date(
    date_breaks = "7 days",
    labels = scales::label_date("%b %d")) +
  geom_histogram(aes(x = doy), binwidth = 1) +
  facet_wrap(vars(year), nrow = 2) +
  labs(
    y = "Surveys",
    x = "Day of year") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/doy.png",
  width = 3.5, height = 3.5, units = "in", dpi = 300)

doy <-
  visits |> 
  mutate(doy = yday(date)) |> 
  select(point_id, visit, year, doy) |> 
  right_join(
    birds |> expand(point_id, year, visit),
    by = c("point_id", "year", "visit")) |> 
  pivot_wider(
    names_prefix = 'doy_',
    names_from = c(visit),
    values_from = doy,
    names_sep = '_') |> 
  unite("point_year", c(point_id, year), remove = T) |> 
  arrange(point_year)

## Combine all obs covariates ----
obs_covs <-
  list(
    observer = observer |> select(-point_year),
    temp = temp |> select(-point_year),
    wind = wind |> select(-point_year),
    precip = precip |>  select(-point_year),
    cloud = cloud |> select(-point_year),
    start_sun = start_sun |> select(-point_year),
    doy = doy |> select(-point_year))

# Generate site covariates ------------------------------------------------

#268 ptyrs
#548 ptyrs

# Import
# NLCD
# Lidar
# Hydrology

site_covs <-
  bind_rows(
    mutate(points, year = 2021),
    mutate(points, year = 2022)) |> 
  # Only visited points
  filter(point_id %in% unique(visits$point_id)) |> 
  # Add plot info
  left_join(plots, by = "plot_name") |> 
  select(point_id, year, plot_type, ownership) |> 
  # Add basal area
  left_join(
    veg |> 
    distinct(point_id, basal_area) |> 
    summarize(
      basal = mean(basal_area, na.rm = T), .by = point_id) |> 
    filter(!is.na(basal)),
  by = "point_id") |> 
# Add provided covariates
  
# Factorize
  mutate(
    plot_type = factor(plot_type), 
    ownership = factor(ownership), 
    year = factor(year)) |> 
# Duplicate for site-year
arrange(point_id, year) 

umfs <- list()

for(i in 1:length(focal_spp)) {
  umfs[i] <-
    unmarkedFramePCount(
    y = 
      abundance |> 
      filter(species == focal_spp[i]) |> 
      select(-species) |> 
      select(-point_year),
    siteCovs = site_covs,
    obsCovs = obs_covs)
}

names(umfs) <- focal_spp

write_rds(
  umfs, 
  "data/processed/umfs.rds")
