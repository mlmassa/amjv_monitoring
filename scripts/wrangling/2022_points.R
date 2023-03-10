# This script maps and generates a kml file of 2022 points to survey
# We used this during the 2022 field season

library(sf)
library(tmap)
library(tidyverse)

# Import data -------------------------------------------------------------

read_rds("data/raw/amjv_data.rds") %>% 
  list2env(.GlobalEnv)

# Import all points
points22 <-
  read_sf("data/raw/GoshenEffectsYr2.shp") %>% 
    select(point_id = Name) %>% 
    st_zm() %>% 
    mutate(plot_name = "Goshen Effects") %>% 
  bind_rows(
    points %>% 
      st_as_sf(
        coords = c("long", "lat"), 
        crs = 4269) %>% 
      filter(!plot_name %in% c("Goshen Effects", "Potts", "TNC Effects")) %>% 
      st_transform(crs = 26917)) %>% 
  arrange(point_id)

# Import blocks
sampling_units <-
  st_read("data/raw/SelectedPSU_Grids.shp") %>% 
  st_transform(crs = 26917) %>% 
  select(
    id = OBJECTID,
    ownership = Own_Type,
    name = Name) %>% 
  mutate(
    name = 
      str_replace_all(name, pattern = "_", replacement = " "),
    ownership = 
      case_when(
        str_detect(ownership, "NGO") ~ "NGO",
        str_detect(ownership, "FED") ~ "Federal",
        str_detect(ownership, "STATE") ~ "State")) %>% 
  filter(
    !name %in% 
      c("State C", "NGO Over", "NGO D", "Fed A", "Fed D", "Fed F"))

# Map points --------------------------------------------------------------

tmap_mode("view")

tm_basemap(
  c("Esri.WorldTopoMap",
    "Esri.WorldImagery")) +
  
  tm_shape(sampling_units) +
  tm_polygons(
    col = "ownership",
    palette = c("#49D4B9", "#D1B7E1", "#FFFF99"),
    alpha = 0.5,
    id = "name",
    popup.vars = 
      c("Plot name" = "name", 
        "Ownership" = "ownership")) +
  
  tm_shape(points22) +
  tm_dots(
    size = 0.05,
    col = "black",
    border.col = "black",
    popup.vars = 
      c("Point ID" = "point_id",
        "Plot name" = "plot_name"),
    clustering = FALSE)

# Write file --------------------------------------------------------------
st_write(
  points22 %>% select(Name = point_id), 
  "data/processed/points_22.kml", 
  driver = "kml",
  delete_dsn = TRUE)
