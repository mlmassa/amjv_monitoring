# Setup -------------------------------------------------------------------

# Load required packages
library(tidyverse)

# Import data
read_rds('data/raw/amjv_data.rds') |> 
  list2env(.GlobalEnv)

# Only retain ARU information
#rm(count_data, metadata, species_list, visits, veg_data)

# Audiomoth status --------------------------------------------------------

# Produce a list of ARUs and most recent action
audiomoth_activity |> 
  group_by(deployment_id) |> 
  # Get most recent action
  filter(
    date == max(date),
    str_detect(deployment_id, 'BATH2022')) |>
  select(deployment_id, action) |> 
  ungroup()

# Get deployment + point
arus_deployed <-
  audiomoth_activity |> 
  distinct(point_id, deployment_id) |> 
  left_join(
    audiomoths |> 
    select(deployment_id, recorder_id, card_id)) |> 
  # Remove undeployed ARUs
  filter(!is.na(point_id)) |> 
  # Add year
  mutate(
    year = if_else(str_detect(deployment_id, 'BATH2021'), 2021, 2022))

# Show errors
duplicates <-
  arus_deployed |> 
  count(deployment_id) |> 
  filter(n!=1) |> 
  pull(deployment_id)

arus_deployed |> 
  filter(deployment_id %in% duplicates) |> 
  arrange(deployment_id)

audiomoth_activity |> 
  filter(deployment_id %in% duplicates) |> 
  arrange(deployment_id)

# Generate ARU report -----------------------------------------------------

aru_report <-
  arus_deployed |> 
  left_join(
    audiomoth_activity |> 
    select(deployment_id, action, date) |> 
    mutate(date = as.character(date)) |> 
    pivot_wider(
      names_from = action, 
      values_from = date)) |> 
    select(
      deployment_id, recorder_id, card_id, point_id,
      year, deployed, collected, sd_removed, downloaded, lost) |> 
  left_join(points) |> 
  left_join(plots)

# Save ARU report
write_csv(aru_report, 'output/amrg_aru_report.csv')

# Veg data ----------------------------------------------------------------
write_csv(veg_data, 'output/amrg_veg_report.csv')

# Aru plot summary --------------------------------------------------------

arus_deployed |> 
  left_join(points) |> 
  left_join(plots) |> 
  group_by(plot_type, year) |> 
  summarize(n = length(unique(point_id)))

# Bear points -------------------------------------------------------------

audiomoth_activity |> 
  # Where did we mention bears
  filter(
    str_detect(aru_note, '[Bb]ear') & !str_detect(aru_note, '[Bb]earing')) |> 
  group_by(year(date)) |> 
  distinct(point_id) |> 
  count()
