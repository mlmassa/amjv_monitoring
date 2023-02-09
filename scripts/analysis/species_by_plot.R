
# setup -------------------------------------------------------------------

# Load required packages

library(leaflet)
library(sf)
library(tmap)
library(lubridate)
library(tidyverse)

# Import data

# Replace with wherever you saved your RDS from assign_event_id.R
read_rds('data/raw/amjv_data.rds') %>% 
  list2env(.GlobalEnv)


# Create untidy dataset ---------------------------------------------------

untidy_data <-
  count_data %>% 
    left_join(visits) %>% 
    left_join(points) %>% 
    left_join(plots)


# Get effort by plot type/block -------------------------------------------

untidy_data %>% 
  group_by(plot_type, year(date)) %>% 
  summarize(
    points_surveyed = length(unique(point_id)),
    species = length(unique(species)))


# Solo counts -------------------------------------------------------------

visits %>% 
  group_by(point_id, year(date)) %>% 
  summarize(visits = length(visit)) %>% 
  filter(visits == 1)


# Histogram of #species per count -----------------------------------------

untidy_data %>% 
  group_by(visit_id) %>% 
  summarize(
    species = length(unique(species))) %>% 
  left_join(visits %>% select(visit_id, point_id)) %>% 
  left_join(points %>% select(point_id, plot_name)) %>% 
  left_join(plots %>% select(plot_name, ownership)) %>% 
  ggplot(aes(species)) +
  geom_histogram(binwidth = 1)

# Boxplot of #species by ownership
untidy_data %>% 
  group_by(visit_id) %>% 
  summarize(
    species = length(unique(species))) %>% 
  left_join(visits %>% select(visit_id, point_id)) %>% 
  left_join(points %>% select(point_id, plot_name)) %>% 
  left_join(plots %>% select(plot_name, ownership)) %>% 
  ggplot(aes(x = ownership, y = species)) +
  geom_boxplot()



# Get effort by plot ------------------------------------------------------

# Generate a tibble of:
# Plot, npoints surveyed, nARUs deployed, nspecies, nWOTH, nCERW, nGWWA


summary_plots <-
  
  # Plots, points, species, birds
  untidy_data %>% 
    group_by(plot_name) %>% 
    summarize(
          points_surveyed = length(unique(point_id)),
          species = length(unique(species))) %>% 
  
  # Add number of focal species per plot
  left_join(
    untidy_data %>% 
      group_by(plot_name, species) %>% 
      summarize(count = n()) %>% 
      filter(species %in% c('WOTH', 'GWWA', 'CERW')) %>% 
      pivot_wider(
        names_from = species,
        values_from = count),
    by = 'plot_name') %>% 
  
  # Add number of ARUs deployed per plot
  left_join(
    audiomoth_activity %>% 
      filter(action == 'deployed') %>% 
      select(point_id) %>% 
      unique() %>% 
      # Add point
      left_join(
        points %>% 
          select(point_id, plot_name),
        by = 'point_id') %>% 
      # Add plot
      left_join(
        plots %>% 
          select(plot_name),
        by = 'plot_name') %>% 
      group_by(plot_name) %>% 
      # Get ARUs per plot
      summarize(
        ARUs = n()),
    by = 'plot_name') %>% 
  
  # Replace NA with 0
  mutate_all(~replace(., is.na(.), 0))


# Get most abundant species -----------------------------------------------

visits %>% filter(visit == 2) %>% select(visit_id) %>% count()
# 118 visits round 1
# 91 visits round 2

abundance <-
  untidy_data %>% 
    group_by(visit_id, visit, species, .drop = FALSE) %>% 
    # Get count of each species on each visit
    summarize(
      count = (n()),
      .groups = 'drop') %>% 
    group_by(species, visit) %>% 
    # Get total abundance, number of counts appearing on
    summarize(
      abundance = sum(count),
      counts = n()) %>% 
    # Summarize presence
    mutate(
      total_counts = if_else(visit == 1, 118, 91),
      prop_counts_present = counts/total_counts,
      avg_per_count = abundance/total_counts,
      avg_if_present = abundance/counts)

# Which species were ONLY detected round 2?
abundance %>% 
  ungroup %>% 
  select(species, visit, prop_counts_present) %>% 
  pivot_wider(names_from = visit, values_from = prop_counts_present, names_prefix = 'visit_') %>% 
  replace_na(list(visit_1 = 0, visit_2 = 0)) %>% 
  filter(visit_1 == 0)

# Get top 15 spp
top_15 <-
  abundance %>% 
  summarize(total_abund = sum(abundance)) %>% 
  arrange(-total_abund) %>% 
  slice_head(n = 20) %>% 
  pull(species)

# Get top species to graph
abundance %>% 
  filter(species %in% top_15) %>% 
  ggplot(
    aes(
      x = reorder(species, -prop_counts_present), 
      y = prop_counts_present, 
      fill = factor(visit))) +
  geom_bar(
    stat = 'identity', 
    width = 0.75,
    color = 'black',
    position = position_dodge(width = 0.75)) +
  labs(
    x = 'Species',
    y = 'Prop. counts present',
    fill = 'Visit') +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1),
    expand = c(0,0)) +
  scale_fill_manual(values = c('grey20', 'grey60')) +
  theme(axis.title = element_text(size = 14),
      plot.title = element_text(size = 16),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.75),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = c(0.9, 0.8),
      panel.background = element_rect(fill = 'white'),
      panel.grid = element_line(color = "#bfbfbf"),
      panel.grid.minor = element_line(linetype = 'dashed'),
      axis.line = element_line(color = 'black'))

  

# Get species by plot -----------------------------------------------------

# Generate a list with a tibble of species for each plot

# Map function:

species_by_plot <-
  map(
    1:length(unique(plots$plot_name)),
    function(i) {
      
      # Species by plot name:
        plot_split <-
          count_data %>% 
            select(visit_id, species) %>% 
            left_join(
              visits %>% 
                select(point_id, visit_id), 
              by = 'visit_id') %>% 
            left_join(
              points %>% 
                select(point_id, plot_name), 
              by = 'point_id') %>% 
          
          # Separate by plot
          filter(plot_name == unique(plots$plot_name)[i]) %>% 
          
          # Add species names
          left_join(
            species_list %>% 
              select(sci_name, common_name, banding_code, taxon_order), 
            by = c('species' = 'banding_code')) %>% 
          
          # Order and generate compact list of unique species
          arrange(taxon_order) %>% 
          select(common_name, sci_name) %>% 
          unique()
    }) %>% 
  # Name list elements by plot
  set_names(unique(plots$plot_name))


# Make an interactive map of points ---------------------------------------


# Convert the stupid datetime

date(visits$start_time) <- date(visits$date)

# Get summary by point

map_data <- 
  untidy_data %>% 
    # Get species count
    group_by(point_id) %>% 
      summarize(
        total_species = length(unique(species)),
        WOTH = length(event_id[species == 'WOTH']),
        CERW = length(event_id[species == 'CERW'])) %>% 
    # Add plot name, lat, long
    right_join(points, by = 'point_id') %>% 
    right_join(plots, by = 'plot_name') %>% 
    # Add ARU number at point
    right_join(
      audiomoth_activity %>% 
        filter(action == 'deployed') %>% 
        select(point_id, recorder_id) %>% 
        unique(),
      by = 'point_id') %>% 
  # Make it spatial
  st_as_sf(
    coords = c('long', 'lat')#,
    #crs = 4326
    )

# Import blocks

sampling_units <-
  st_read('data/raw/SelectedPSU_Grids.shp') %>% 
  #st_transform(crs = 4326) %>% 
  select(id = OBJECTID,
         ownership = Own_Type,
         name = Name) %>% 
  mutate(
    name = 
      str_replace_all(name, pattern = '_', replacement = ' '),
    ownership = 
      case_when(
        str_detect(ownership, 'NGO') ~ 'NGO',
        str_detect(ownership, 'FED') ~ 'Federal',
        str_detect(ownership, 'STATE') ~ 'State')) %>% 
  filter(
    !name %in% 
      c('State C', 'NGO Over', 'NGO D', 'Fed A', 'Fed D', 'Fed F')) %>% 
  left_join(
    summary_plots %>% 
      select(plot_name, points_surveyed, species, ARUs),
    by = c('name' = 'plot_name')
  )


# Map of points, colored by species richness

tmap_mode('view')
tm_view(
  set.view = c(-79.799709, 37.994635, 11)) +

tm_basemap(
  c('Esri.WorldTopoMap',
    'Esri.WorldImagery')) +
  
  tm_shape(sampling_units) +
  tm_polygons(
    col = 'ownership',
    palette = c('#49D4B9', '#D1B7E1', '#FFFF99'),
    alpha = 0.5,
    id = 'name',
    popup.vars = 
      c('Plot name' = 'name', 
        'Ownership' = 'ownership', 
        'Points surveyed' = 'points_surveyed',
        'Total species' = 'species',
        'ARUs deployed' = 'ARUs')) +
  
  tm_shape(map_data) +
  tm_dots(
    size = 0.05,
    col = 'total_species',
    border.col = 'black',
    popup.vars = 
      c('Plot name' = 'plot_name',
        'Ownership' = 'ownership',
        'Monitoring type' = 'plot_type',
        'Species' = 'total_species', 
        'ARU ID' = 'recorder_id'),
    clustering = FALSE)
  

# Map WOTH by point -------------------------------------------------------

# Map of points, colored by WOTH detections (0-4)

tmap_mode('view')
tmap_options(check.and.fix = T)

tm_basemap(
  c('Esri.WorldTopoMap',
    'Esri.WorldImagery')) +
  
  tm_shape(map_data) +
  tm_dots(
    size = 0.1,
    col = 'WOTH',
    palette = 'BuGn',
    popup.vars = c('plot_name',
                   'plot_type',
                   'total_species',
                   'WOTH',
                   'CERW', 
                   'recorder_id'),
    title = 'WOTH detections',
    colorNA = 'grey40',
    textNA = 'Not surveyed',
    clustering = FALSE) 

  
# Most abundant species ---------------------------------------------------

untidy_data %>% 
  group_by(visit, species) %>% 
  summarize(
    count = length(unique(event_id))) %>% 
  pivot_wider(names_from = visit, values_from = count) %>% 
  arrange(-`1`)
    