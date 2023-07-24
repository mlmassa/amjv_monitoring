# Format data for occupancy analysis

# Setup -------------------------------------------------------------------

library(tidyverse)
library(suncalc)
library(unmarked)

theme_set(
  theme_bw(base_size = 11) +
  theme(
    axis.ticks = element_line(color = "gray"),
    panel.border = element_rect(fill = NA, color = "gray"),
    panel.grid.major = element_line(color = "gray95"),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()))

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

spat_covs <-
  read_rds("data/processed/spat_covs.rds") |> 
  sf::st_drop_geometry()

# Fix the time and get sunrise
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
  filter(
    flyover == 0,
    dist_interval < 5) |> 
  # Only birds within count circle
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
  facet_wrap(vars(species), scales = "free_y") +
  labs(x = "Abundance (birds/survey)", y = "Surveys") +
  coord_cartesian(expand = T, clip = "off")

write_rds(detections, "data/processed/detections.rds")

ggsave(
  "output/plots/focal_abund.png",
  width = 6, height = 5, units = "in", dpi = "retina")

# Plot occupancy
detections |> 
  filter(!is.na(abund)) |> 
  ggplot() + 
  geom_bar(aes(x = factor(det))) + 
  facet_wrap(vars(species)) +
  labs(x = "Detected", y = "Surveys") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/focal_det.png",
  width = 6, height = 5, units = "in", dpi = "retina")

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

write_rds(occupancy, "data/processed/occupancies.rds")

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
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

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
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

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
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

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
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

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
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

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
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

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
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

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

# Import site covariates
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
# Add spatial covariates (lsm, area, stream, forest...)
  left_join(
    sf::st_drop_geometry(spat_covs),
    by = "point_id") |> 
# Factorize
  mutate(
    point_id = factor(point_id),
    plot_type = factor(plot_type), 
    plot_name = factor(plot_name),
    ownership = factor(ownership), 
    year = factor(year),
    d_stream = units::drop_units(d_stream),
    ecozone = factor(ecozone)) |> 
# Duplicate for site-year
arrange(point_id, year) 

# Plot site covariates

# Plot type
slice_head(site_covs, n = 1, by = point_id) |> 
  ggplot() +
  geom_bar(aes(x = plot_type)) +
  labs(x = "Plot type", y = "Points") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/plottype.png",
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

# Ownership
slice_head(site_covs, n = 1, by = point_id) |> 
  ggplot() +
  geom_bar(aes(x = ownership)) +
  labs(x = "Ownership", y = "Points") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/ownership.png",
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

# Basal area
slice_head(site_covs, n = 1, by = point_id) |> 
  ggplot() +
  geom_histogram(aes(x = basal), binwidth = 10) +
  labs(x = "Basal area", y = "Points") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/basal.png",
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

# Ecozone
slice_head(site_covs, n = 1, by = point_id) |> 
  ggplot() +
  geom_bar(aes(x = forcats::fct_infreq(ecozone))) +
  labs(x = "Ecozone", y = "Points") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_flip(clip = "off")

ggsave(
  "output/plots/ecozone.png",
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

# Forest buffers
slice_head(site_covs, n = 1, by = point_id) |> 
  select(point_id, for_1500, for_5000) |> 
  pivot_longer(cols = c(for_1500, for_5000)) |> 
  mutate(
    name = if_else(name == "for_1500", "1.5 km", "5 km")) |> 
  ggplot() +
  geom_histogram(aes(x = value), binwidth = .01) +
  facet_wrap(~name, nrow = 2) +
  labs(x = "Prop. forest within buffer", y = "Points") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/forest.png",
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

# Core patch area
slice_head(site_covs, n = 1, by = point_id) |> 
  ggplot() +
  geom_histogram(aes(x = core), bins = 10) +
  labs(x = "Core patch area (ha)", y = "Points") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/core.png",
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

# Patch perim-area ratio
slice_head(site_covs, n = 1, by = point_id) |> 
  ggplot() +
  geom_histogram(aes(x = para), bins = 20) +
  labs(x = "Patch perimeter:area ratio", y = "Points") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/para.png",
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

# Distance to stream
slice_head(site_covs, n = 1, by = point_id) |> 
  ggplot() +
  geom_histogram(aes(x = d_stream), binwidth = 10) +
  labs(x = "Distance to stream (m)", y = "Points") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/dstream.png",
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

# Edge density 1 mile
slice_head(site_covs, n = 1, by = point_id) |> 
  ggplot() +
  geom_histogram(aes(x = ed_buff), binwidth = 2) +
  labs(x = "Forest edge density (m/ha)", y = "Points") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/ed.png",
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

# Mean canopy height within 50m
slice_head(site_covs, n = 1, by = point_id) |> 
  ggplot() +
  geom_histogram(aes(x = canopy_h), binwidth = 1) +
  labs(x = "Mean canopy height (m) within 50 m", y = "Points") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/canopy.png",
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

# Heterogeneity (CV of canopy height within 100m)
slice_head(site_covs, n = 1, by = point_id) |> 
  ggplot() +
  geom_histogram(aes(x = canopy_cv), binwidth = 0.05) +
  labs(x = "CV of canopy height (m) within 100 m", y = "Points") +
  coord_cartesian(expand = T, clip = "off")

ggsave(
  "output/plots/canopycv.png",
  width = 3.5, height = 3.5, units = "in", dpi = "retina")

# Make UMFs ---------------------------------------------------------------

umfs <- list()

for(i in 1:length(focal_spp)) {
  umfs[[i]] <-
    unmarkedFramePCount(
    y = 
      abundance |> 
      filter(species == focal_spp[i]) |> 
      select(-species) |> 
      select(-point_year),
    siteCovs = site_covs |> sf::st_drop_geometry(),
    obsCovs = obs_covs |> sf::st_drop_geometry())
}

names(umfs) <- focal_spp

write_rds(
  umfs, 
  "data/processed/umfs.rds")
