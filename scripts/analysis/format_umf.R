# Format data for occupancy analysis

# Setup -------------------------------------------------------------------

library(tidyverse)
library(unmarked)

# Import data
counts <- read_rds("data/raw/amjv_data.rds")$count_data
visits <- 
  read_rds("data/raw/amjv_data.rds")$visits %>% 
  mutate(
    year = lubridate::year(date),
    # Create combined pointxyear variable (to serve as site for single-season model)
    point_year = str_c(point_id, "_", 1 + year - min(year)))
points <- read_rds("data/raw/amjv_data.rds")$points
plots <- read_rds("data/raw/amjv_data.rds")$plots

veg <- read_rds("data/raw/amjv_data.rds")$veg_data

# Generate untidy dataset:
birds <-
  counts %>% 
  left_join(visits) %>% 
  left_join(points) %>% 
  left_join(plots)

# Set parameters
focal_spp <- c("WOTH", "CERW")

# Generate detection table ------------------------------------------------

# Create annual point visit table
annual_points <-
  visits %>% 
  # Get unique points for each year
  group_by(year) %>%
  select(point_id, year) %>% 
  distinct() %>% 
  # Mark these as visited
  mutate(visited = 1) %>% 
  ungroup() %>% 
  # Get full table of points per year
  complete(year, nesting(point_id), fill = list(visited = 0)) %>% 
  arrange(point_id, year)

# Generate detection history
# Problematic, needs work
detections <-
  birds %>% 
  # No flyover detections: birds in habitat only
  filter(flyover == 0) %>% 
  # For each species at each pointxyearxvisit...
  group_by(point_id, year, point_year, visit, species) %>% 
  # Get abundance
  summarize(abund = n(), .groups = "drop") %>% 
  # Add empty combinations not present
  right_join(
    birds %>% expand(point_id, year, visit, species),
    by = c("point_id", "year", "visit", "species")) %>% 
  # Add whether pointxyear combo was visited
  left_join(annual_points) %>%
  # Add detection (1/0 if visited, keep NA if not)
  mutate(det = 
    case_when(
      is.na(abund) & visited == 1 ~ 0,
      abund > 0 ~ 1)) %>% 
  # Retain only focal species
  filter(species %in% focal_spp) %>% 
  # Cut irrelevant variables
  select(point_year, point_id, year, visit, species, det) %>% 
  # Make wide data with column for each year-visit detection
  # detection is "y"
  pivot_wider(
    names_prefix = 'y_', 
    names_from = c(visit),
    values_from = det,
    names_sep = '_') %>% 
  arrange(point_id, year) 

# Generate observation covariates -----------------------------------------

# Format:
# siteCovs: single (n_sitexyear) row x n_covariates column DF 
# obsCovs: list of (n_sitexyear) row x n_visits column DFs

visits %>% 
  select(point_id, visit, year, observer) %>% 
  mutate(observer = as.factor(observer)) %>% 
  pivot_wider(
    names_prefix = 'observer_',
    names_from = c(visit),
    values_from = observer,
    names_sep = '_') %>% 
# Standardize formatting
arrange(point_id, year) %>% 
mutate(year = year - 2013) %>% 
unite('grts_year', c(grts, year)) %>% 
full_join(
  detections %>% select(grts_year) %>% unique()) %>% 
column_to_rownames(var = 'grts_year')



