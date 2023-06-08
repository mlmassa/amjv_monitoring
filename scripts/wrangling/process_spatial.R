# Process and reprocess spatial covariates
# Also make a map of the study area

# Setup -------------------------------------------------------------------

library(tidyverse) # Data manipulation, graphing, etc.
library(sf) # Handle spatial data
library(tigris) # Fetch county/state/road/water shapefiles
options(tigris_use_cache = TRUE) # Don't re-download shapefiles
library(terra) # Handle raster files
library(FedData) # Get NLCD
library(landscapemetrics) # Calculate landscape metrics

# Data import -------------------------------------------------------------

# Import points and blocks
points <-
  # Import point locations
  read_rds("data/raw/amjv_data.rds")$points |> 
  # Keep most things as longlat unprojected (EPSG 4326)
  st_as_sf(coords = c("long", "lat"), crs = 4326) |> 
  # Keep only visited points
  filter(point_id %in% read_rds("data/raw/amjv_data.rds")$visits$point_id)

plots <-
  # Import plot data ("blocks")
  read_rds("data/raw/amjv_data.rds")$plots |> 
  # Add spatial outline for plot
  left_join(
    st_read("data/raw/SelectedPSU_Grids.shp") |> 
      # Fix plot names
      mutate(plot_name = str_replace_all(Name, "_", " ")) |> 
      # Drop unneeded information
      select(plot_name, geometry) |> 
      st_transform(crs = st_crs(points))) |> 
  # Remove unvisited plots
  filter(plot_name %in% points$plot_name) |> 
  st_as_sf()

# Get large bounding box (all points)
bbox_large <-
  points |> 
  left_join(st_drop_geometry(plots)) |> 
  # Get 6km buffer
  st_buffer(dist = 6000) |> 
  st_bbox()

# Get small bbox (no Potts)
bbox_small <-
  points |> 
  filter(plot_name != "Potts") |> 
  left_join(st_drop_geometry(plots)) |> 
  # Get 3km buffer
  st_buffer(dist = 3000) |> 
  st_bbox()

# Ecozone -----------------------------------------------------------------

ecozone_xwalk <-
  readxl::read_xls("data/raw/zone_codes.xls")

ecozone_rast <-
  rast("data/raw/Cropped_Raster/AppRidges_cropped.tif") |> 
  classify(rcl = as.matrix(select(ecozone_xwalk, Value, New_Value)))

names(ecozone_rast) <- "New_Value"

# Add crosswalk level name to raster value
levels(ecozone_rast) <- 
  distinct(ecozone_xwalk, New_Value, GWNF_model) |> 
  arrange(New_Value)

# Re-run extraction
ecozone <-
  points |> 
  mutate(
    extract(
      x = ecozone_rast,
      y = points |> 
        st_transform(crs = st_crs(ecozone_rast)) |> 
        vect(),
      ID = F))

# Forest cover ------------------------------------------------------------

# # Download NLCD (only need to run once)
# get_nlcd(
#   template = st_as_sfc(bbox_large),
#   label = "amjv",
#   year = 2019,
#   dataset = "landcover",
#   extraction.dir = "data/raw/",
#   raster.options = c(overwrite = TRUE))

# Import NLCD
nlcd <- 
  rast("data/raw/amjv_NLCD_Land_Cover_2019.tif")

# Create forest reclassification matrix
# Lumping all forest (41 Deciduous, 42 Evergreen, 43 Mixed)
nlcd_forest <-
  classify(
    nlcd,
    tibble(from = sort(unique(values(nlcd)))) %>% 
      mutate(to = if_else(from %in% c(41:43), 1, 0)) %>% 
      as.matrix())

# Transform and buffer points for extraction
radius <-
  1609 # 1 mi
  #2000 # 2 km
  #5000 # 5 km
  
points_buff <- 
  st_transform(points, crs = st_crs(nlcd)) |> 
  st_buffer(dist = radius) 

# Perform extraction
forest <-
  points |> 
  mutate(
    extract(
      x = nlcd_forest,
      y = vect(points_buff),
      fun = mean,
      ID = F)) |> 
  rename(for_1mi = Class)

# Landscape metrics -------------------------------------------------------

# Core area: >200m from edge

# Check if valid for landscapemetrics
check_landscape(nlcd_forest)

# Extract patch metrics
metrics <-
  extract_lsm(
    landscape = nlcd_forest,
    y = points |> 
      st_transform(crs = st_crs(nlcd_forest)) |> 
      select(geometry),
    extract_id = points$point_id,
    what = c(
      "lsm_p_core", # Core patch area
      "lsm_p_para" # Perimeter-area ratio of patch
    ),
    #200m/30m resolution = 7 cells
    edge_depth = 7) |> 
  pivot_wider(
    id_cols = extract_id,
    names_from = metric,
    values_from = value) |> 
  rename(point_id = extract_id)

# Point-buffer metrics
# Get edge density within buffer radius
ed <-
  map(
    # For each buffered point...
    1:nrow(points_buff),
    # Get the forest edge density within that buffer
    ~lsm_l_ed(
        mask(nlcd_forest, vect(points_buff[.x, ])) |> 
        crop(vect(points_buff[.x, ])),
      count_boundary = F)) |> 
  bind_rows() |> 
  mutate(point_id = points_buff$point_id) |> 
  select(point_id, ed_1mi = value)

# Distance to stream ------------------------------------------------------

# Stream linear data from VA DCR
streams <-
  st_read("data/raw/Rivers_(2022_Final_WQA_IR_Assessment).shp") |> 
  st_transform(crs = st_crs(points)) |> 
  st_crop(st_as_sfc(bbox_large))

# Water area features
water_area <-
  bind_rows(
    # Focal areas
    area_water(
      state = "VA", 
      county = "Rockbridge"),
    area_water(
      state = "VA", 
      county = "Bath"),
    area_water(
      state = "VA", 
      county = "Alleghany"),
    # Other areas that will appear on map
    area_water(
      state = "VA", 
      county = "Covington"),
    area_water(
      state = "VA", 
      county = "Botetourt"),
    area_water(
      state = "VA", 
      county = "Craig"),
    area_water(
      state = "WV", 
      county = "Greenbrier"),
    area_water(
      state = "WV", 
      county = "Monroe"),
    area_water(
      state = "WV", 
      county = "Pocahontas")) |>
  st_transform(crs = st_crs(points))

# Calculate distance to stream
dist_stream <-
  st_distance(points, streams)

rownames(dist_stream) <- points$point_id

# Calculate distance to water body
dist_area <-
  st_distance(points, water_area)

rownames(dist_area) <- points$point_id

d_stream <-
  bind_rows(
    as_tibble(t(dist_stream)),
    as_tibble(t(dist_area))) |> 
  # Get minimum distance to stream OR water body
  summarize(
    across(1:134, ~min(.x, na.rm = T))) |> 
  pivot_longer(
    cols = 1:134,
    names_to = "point_id",
    values_to = "d_stream")

# Canopy height -----------------------------------------------------------

canopy_rast <-
  rast("data/raw/canopyheight.tif")

# Canopy height (avg. within 50m)
canopy <-
  points |> 
  mutate(
    extract(
      x = canopy_rast,
      y = 
        points |> 
        st_transform(crs = st_crs(canopy_rast)) |> 
        st_buffer(dist = 50) |> 
        vect(),
      fun = mean,
      ID = F)) 

# Export spatial covariates -----------------------------------------------

spat_covs <-
  st_drop_geometry(ecozone) |> 
  left_join(st_drop_geometry(forest)) |> 
  left_join(metrics) |> 
  left_join(d_stream) |> 
  left_join(ed) |> 
  left_join(canopy) |> 
  rename(
    ecozone = GWNF_model,
    h_canopy = canopyheight)

write_rds(spat_covs, "data/processed/spat_covs.rds")

# Extra map data ----------------------------------------------------------

# Add additional packages
library(elevatr) # Fetch elevation data
library(tidyterra)  # Make rasters work with ggplot

## State borders ----
states <-
  states() |> 
  filter(NAME %in% c("Virginia", "West Virginia")) |> 
  st_transform(crs = st_crs(points))

## County borders ----
counties <-
  counties(state = states$STUSPS) |> 
  st_transform(crs = st_crs(points))

## Main roads ----
roads <-
  bind_rows(
    primary_secondary_roads(state = "VA"),
    primary_secondary_roads(state = "WV")) |> 
  st_transform(crs = st_crs(points)) |> 
  st_crop(bbox_large)

## Water (linear) ----
water_line <-
  bind_rows(
    # Focal areas
    linear_water(
      state = "VA", 
      county = "Rockbridge"),
    linear_water(
      state = "VA", 
      county = "Bath"),
    linear_water(
      state = "VA", 
      county = "Alleghany"),
    # Other areas that will appear on map
    linear_water(
      state = "VA", 
      county = "Covington"),
    linear_water(
      state = "VA", 
      county = "Botetourt"),
    linear_water(
      state = "VA", 
      county = "Craig"),
    linear_water(
      state = "WV", 
      county = "Greenbrier"),
    linear_water(
      state = "WV", 
      county = "Monroe"),
    linear_water(
      state = "WV", 
      county = "Pocahontas")) |>
  st_transform(crs = st_crs(points))

## Hillshade ----
# Set hillshade palette
grays <- 
  colorRampPalette(c("gray65", "gray98"))

# Get vector of hillshade colors
pal_gray <- 
  grays(1000)

# Download elevation raster
elev_large <-
  get_elev_raster(
    as_Spatial(st_as_sfc(bbox_large)), 
    neg_to_na = T, # Large negative numbers: eliminate
    z = 11) |> # Zoom level
  # Downloads as raster; convert to terra
  rast()

# Create hillshade effect
hill_large <- 
  shade(
    terrain(elev_large, "slope", unit = "radians"), 
    terrain(elev_large, "aspect", unit = "radians"), 
    30, 270)

# Scale raster to # of colors in gray palette
index_large <- 
  hill_large |> 
  mutate(
    index_col = scales::rescale(hillshade, to = c(1, length(pal_gray)))) |>
  mutate(index_col = round(index_col)) |> 
  pull(index_col)

# Set color palette of hillshade
vector_cols_large <- 
  pal_gray[index_large]

# Plot --------------------------------------------------------------------

# Which to plot
bbox <- bbox_small
bbox <- bbox_large

# Set water color
water_col <- "cornflowerblue"

# Create plot labels
plot_labels <-
  points |> 
  group_by(plot_name) |> 
  summarise(geometry = st_union(geometry))

# Attempt to fix label text size.
# Not sure if this works
text_size_pt <- 6
text_size <- text_size_pt * 0.352777777777777
update_geom_defaults("label", list(size = text_size))
update_geom_defaults("text", list(size = text_size))

ggplot() +
  ## Hillshade ----
  geom_spatraster(
    data = hill_large, 
    fill = vector_cols_large, 
    maxcell = Inf,
    alpha = 1) +
  ## Water area ----
  geom_sf(
    data = water_area |> st_transform(crs = st_crs(hill_large)),
    fill = water_col,
    color = water_col,
    alpha = 0.9) +
  ## Water linear ----
  geom_sf(
    data = water_line |> st_transform(crs = st_crs(hill_large)),
    color = water_col,
    linewidth = 0.1,
    alpha = 0.9) +
  geom_sf(
    data = streams |> st_transform(crs = st_crs(hill_large)),
    color = water_col,
    linewidth = 0.1,
    alpha = 0.65) +
  ## Roads ----
  geom_sf(
    data = roads |> st_transform(crs = st_crs(hill_large)),
    # Scale line width width by road type
    aes(linewidth = MTFCC),
    color = "white") +
  scale_linewidth_manual(values = c(0.6, 0.3), guide = "none") +
  ## Counties ----
  geom_sf(
    data = counties |> st_transform(crs = st_crs(hill_large)), 
    fill = NA,
    color = alpha("gray30", 0.8),
    linetype = "dashed",
    linewidth = 0.2) +
  ## US state borders ----
  geom_sf(
    data = states |> st_transform(crs = st_crs(hill_large)), 
    fill = NA,
    linewidth = 0.25,
    color = "gray30") +
  ## Plots ----
  geom_sf(
    data = plots |> st_transform(crs = st_crs(hill_large)), 
    aes(color = ownership),
    linewidth = 0.25,
    fill = NA) +
  ## Plot labels ----
  ggrepel::geom_label_repel(
    data = plot_labels |> st_transform(crs = st_crs(hill_large)) 
    # |> 
    #   filter(plot_name != "Potts")
    ,
    aes(label = plot_name, geometry = geometry),
    stat = "sf_coordinates",
    size = 2,
    label.padding = 0.2,
    force_pull = 0.85,
    min.segment.length = unit(0.2, "inches")) +
  ## Points ----
  geom_sf(
    data = points |> 
      left_join(st_drop_geometry(plots)) |> 
      st_transform(crs = st_crs(hill_large)),
    aes(shape = plot_type),
    color = "black", 
    size = 0.7) +
  scale_shape_manual(values = c(3, 20)) +
  ## Theme ----
  theme_void() +
  coord_sf(
    xlim = c(bbox$xmin, bbox$xmax),
    ylim = c(bbox$ymin, bbox$ymax),
    expand = F) +
  labs(
    color = "",
    fill = "",
    shape = "",
    linewidth = "") +
  theme(
    panel.background = element_rect(fill = "#F1F3F4"),
    legend.position = "top",
    panel.border = element_rect(fill = NA))

## Save plot ----
ggsave(
  "output/plots/map.png",
  width = 7.5,
  height = 7,
  units = "in",
  dpi = "retina")
