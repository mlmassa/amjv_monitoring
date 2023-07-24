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
  # Get buffer
  st_buffer(dist = 10000) |> 
  st_bbox()

# Get small bbox (no Potts)
bbox_small <-
  points |> 
  filter(plot_name != "Potts") |> 
  left_join(st_drop_geometry(plots)) |> 
  # Get 3km buffer
  st_buffer(dist = 3000) |> 
  st_bbox()

# Get Potts only bbox
bbox_potts <-
  points |> 
  filter(plot_name == "Potts") |> 
  left_join(st_drop_geometry(plots)) |> 
  # Get buffer
  st_buffer(dist = 2000) |> 
  st_bbox()

# Get Gathright bbox
bbox_gathright <-
  points |> 
  left_join(st_drop_geometry(plots)) |> 
  filter(block == "Gathright", plot_name != "Potts") |> 
  # Get buffer
  st_buffer(dist = 3000) |> 
  st_bbox()

# Get Goshen bbox
bbox_goshen <-
  points |> 
  left_join(st_drop_geometry(plots)) |> 
  filter(block == "Goshen", plot_name != "Potts") |> 
  # Get buffer
  st_buffer(dist = 2000) |> 
  st_bbox()

# Ecozone -----------------------------------------------------------------

ecozone_xwalk <-
  readxl::read_xls("data/raw/zone_codes.xls")

ecozone_rast <-
  rast("data/raw/Cropped_Raster/AppRidges_cropped.tif") |> 
  # Reclassify ecozone to fewer categories using crosswalk
  classify(rcl = as.matrix(select(ecozone_xwalk, Value, New_Value)))

names(ecozone_rast) <- "New_Value"

# Add crosswalk level name to raster value
levels(ecozone_rast) <- 
  distinct(ecozone_xwalk, New_Value, GWNF_model) |> 
  arrange(New_Value)

# Run extraction
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
radii <-
  c(1500, # A little less than 1 mile
    5000) # 5 km

# Perform extraction
forest <-
  map(
    radii,
    ~ points |> 
      mutate(
        extract(
          x = nlcd_forest,
          y = points |> 
            st_transform(crs = st_crs(nlcd_forest)) |> 
            st_buffer(dist = .x) |> 
            vect(),
          fun = mean,
          ID = F)) |> 
      mutate(radius = .x) |> 
      select(point_id, forest = Class, radius)) |> 
  list_rbind() |> 
  tibble() |> 
  select(-geometry) |> 
  pivot_wider(
    names_from = radius,
    names_prefix = "for_",
    values_from = forest)

# Landscape metrics -------------------------------------------------------

# Core area: >200m from edge

# Check if valid for landscapemetrics
check_landscape(nlcd_forest)

# Extract patch metrics
metrics <-
  extract_lsm(
    landscape = nlcd_forest,
    y = 
      points |> 
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

# Set edge density radius
points_buff <-
  st_transform(points, crs = st_crs(nlcd)) |>
  st_buffer(dist = 1500)

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
  select(point_id, ed_buff = value)

# Distance to stream ------------------------------------------------------

# Stream linear data from VA DCR
streams <-
  st_read("data/raw/Rivers_(2022_Final_WQA_IR_Assessment).shp") |> 
  st_transform(crs = st_crs(points)) |> 
  st_crop(st_as_sfc(bbox_large))

# Water area features from tigris
# I could use VA DCR but only Moomaw is needed, plus I will use these
# for making a map later
water_area <-
  bind_rows(
    # Focal areas
    area_water(
      state = "VA", 
      county = c("Rockbridge", "Bath", "Alleghany", "Covington", "Botetourt", "Craig")),
    area_water(
      state = "WV", 
      county = c("Greenbrier", "Monroe", "Pocahontas"))) |>
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

canopy_rast_potts <-
  rast("data/raw/canopyheight_potts.tif")

# Canopy height (avg. within 50m)
canopy <-
  bind_rows(
    # Non-Potts points
    filter(points, plot_name != "Potts") |> 
    mutate(
      extract(
        x = canopy_rast,
        y = 
          filter(points, plot_name != "Potts") |> 
          st_transform(crs = st_crs(canopy_rast)) |> 
          st_buffer(dist = 50) |> 
          vect(),
        fun = mean,
        ID = F)),
    # Potts points
    filter(points, plot_name == "Potts") |> 
    mutate(
      extract(
        x = canopy_rast_potts,
        y = 
          points |> 
          filter(plot_name == "Potts") |> 
          st_transform(crs = st_crs(canopy_rast_potts)) |> 
          st_buffer(dist = 50) |> 
          vect(),
        fun = mean,
        ID = F)) |> 
    rename(canopyheight = canopyheight_potts)
  )

cv <-
  function(x, na.rm = T) {sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)}

p_potts <-
  filter(points, plot_name == "Potts")

p_notpotts <-
  filter(points, plot_name != "Potts")

cv_rast <-
  function(points, layer) {
    points |> 
    mutate(canopycv = 
      map(
        1:nrow(points),
        ~ crop(
            layer,
            points[.x ,] |> 
              st_transform(crs = st_crs(layer)) |> 
              st_buffer(dist = 100) |> 
              vect(),
            mask = T) |> 
          values() |> 
          c() |> 
          cv()) |> 
        list_c())
  }

canopy_cv <-
  bind_rows(
    cv_rast(p_potts, canopy_rast_potts),
    cv_rast(p_notpotts, canopy_rast))

canopy <-
  canopy |> 
  left_join(st_drop_geometry(canopy_cv))


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
    canopy_h = canopyheight,
    canopy_cv = canopycv) |> 
  st_drop_geometry() |> 
  select(-geometry)

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
    roads(
      state = "VA", 
      county = c("Rockbridge", "Bath", "Alleghany", "Covington", "Botetourt", "Craig")),
    roads(
      state = "WV",
      county = c("Greenbrier", "Monroe", "Pocahontas"))) |> 
  filter(MTFCC %in% c("S1100", "S1200", "S1400", "S1500", "S1740")) |>
  mutate(
    MTFCC = factor(
      MTFCC, 
      labels = c("Primary", "Secondary", "Local", "4WD", "Logging"))) |> 
  st_transform(crs = st_crs(points)) |> 
  st_crop(bbox_large)

## Water (linear) ----
water_line <-
  bind_rows(
    # Focal areas
    linear_water(
      state = "VA", 
      county = c("Rockbridge", "Bath", "Alleghany", "Covington", "Botetourt", "Craig")),
    linear_water(
      state = "WV", 
      county = c("Greenbrier", "Monroe", "Pocahontas"))) |>
  st_transform(crs = st_crs(points))

## Hillshade ----
# Get vector of hillshade colors
pal_gray <- 
  hcl.colors(1000, "Grays")

# Download elevation raster
elev_large <-
  get_elev_raster(
    as_Spatial(st_as_sfc(bbox_large)), 
    neg_to_na = T, # Large negative numbers: eliminate
    z = 13) |> # Zoom level
  # Downloads as raster; convert to terra
  rast()

# Create hillshade effect
hill_large <- 
  shade(
    terrain(elev_large, "slope", unit = "radians"), 
    terrain(elev_large, "aspect", unit = "radians"), 
    30, 270)

## Protected areas ----
prot_areas <-
  st_read("data/processed/prot_areas.shp") |> 
  st_transform(crs = st_crs(hill_large)) |> 
  filter(
    str_detect(
      loc_nam, 
      "Wildlife Management Area|National Forest|State Forest")) |> 
  mutate(loc_nam = str_remove_all(loc_nam, "T.M. ")) |> 
  st_crop(bbox_large) |> 
  group_by(loc_nam) |> 
  summarize() |> 
  mutate(
    category = case_when(
      str_detect(loc_nam, "Wildlife Management Area") ~ "WMA",
      str_detect(loc_nam, "State Forest") ~ "State Forest",
      str_detect(loc_nam, "National Forest") ~ "National Forest")) |> 
  # These overlap with a NF or are not relevant
  filter(
    category != "State Forest",
    !str_detect(loc_nam, "Rimel|Neola"))

# Plot --------------------------------------------------------------------

# Which to plot
bbox <- 
  counties |> 
  filter(
    NAME %in% c("Highland", "Bath", "Alleghany", "Craig", 
                "Botetourt", "Rockbridge", "Augusta")) |> 
  st_bbox()

bbox <- 
  bbox_gathright

bbox <- 
  bbox_goshen

bbox <- 
  bbox_potts

# Set water color
water_col <- "#0978AB"

# Create plot labels
plot_labels <-
  points |> 
  group_by(plot_name) |> 
  summarise(geometry = st_union(geometry))

# Set color scale of raster, and make it smaller
hill_small <-
  hill_large |> 
  crop(
    bbox |> 
      st_as_sfc() |> 
      vect(),
    mask = T)

# Scale raster to # of colors in gray palette
index_small <- 
  hill_small |> 
  mutate(
    index_col = scales::rescale(hillshade, to = c(1, length(pal_gray)))) |>
  mutate(index_col = round(index_col)) |> 
  pull(index_col)

# Set color palette of hillshade
vector_cols_small <- 
  pal_gray[index_small]

# Set road scales
road_scales <- c(0.1, 0.65, 0.2, 0.2, 0)

gc();ggplot() +
## Hillshade ----
  geom_spatraster(
    data = hill_small,
    fill = vector_cols_small,
    maxcell = Inf,
    alpha = 1) +
## Elevation ----
  geom_spatraster(
    data = elev_large) +
  scale_fill_whitebox_c(
    palette = "gn_yl",
    alpha = 0.5) +
## NF/WMA ----
  # geom_sf(
  #   data = prot_areas |> filter(category == "National Forest"),
  #   fill = alpha("purple", 0.2),
  #   color = "purple",
  #   linewidth = 0.5) +
  # geom_sf(
  #   data = prot_areas |> filter(category == "WMA"),
  #   fill = alpha("yellow", 0.2),
  #   color = "yellow",
  #   linewidth = 0.5) +
## Water area ----
  geom_sf(
    data = water_area,
    fill = water_col,
    color = water_col,
    alpha = 0.9) +
## Water linear ----
  # geom_sf(
  #   data = water_line,
  #   color = water_col,
  #   linewidth = 0.1,
  #   alpha = 0.9) +
  geom_sf(
    data = streams,
    color = water_col,
    linewidth = 0.25,
    alpha = 0.85) +
## Roads ----
  geom_sf(
    data = roads, #|> 
      #filter(MTFCC %in% c("Primary", "Secondary", "Local"), RTTYP != "M"),
    # Scale line width width by road type
    aes(linewidth = MTFCC),
    color = "white") +
  scale_linewidth_manual(values = road_scales, guide = "none") +
## Counties ----
  geom_sf(
    data = counties, 
    fill = NA,
    color = alpha("gray30", 0.8),
    linetype = "dashed",
    linewidth = 0.3) +
## US state borders ----
  geom_sf(
    data = states, 
    fill = NA,
    linewidth = 0.6,
    color = "gray30") +
## Plots ----
  geom_sf(
    data = plots,
    linewidth = 0.5,
    alpha = 0.9,
    color = "yellow",
    fill = NA) +
## Points ----
  geom_sf(
    data = points |>
      left_join(st_drop_geometry(plots)),
    fill = "black",
    color = "white",
    aes(shape = plot_type),
    size = 1.5,
    stroke = 0.25) +
  scale_shape_manual(values = c(24, 21)) +
## Inset maps ----
  # geom_sf(
  #   data = st_as_sfc(bbox_potts),
  #   fill = NA,
  #   col = "orange",
  #   linewidth = 1) +
  # geom_sf(
  #   data = st_as_sfc(bbox_gathright),
  #   fill = NA,
  #   col = "orange",
  #   linewidth = 1) +
  # geom_sf(
  #   data = st_as_sfc(bbox_goshen),
  #   fill = NA,
  #   col = "orange",
  #   linewidth = 1) +
## Scale bar ----
  ggspatial::annotation_scale(
    style = "ticks",
    location = "br") +
  ggspatial::annotation_north_arrow(
    style = ggspatial::north_arrow_orienteering(
      text_col = NA,
      fill = c("black", "black")),
    width = unit(0.3, "cm"),
    height = unit(0.5, "cm"),
    location = "tr") +
## Theme ----
  theme_void() +
  guides(fill = "none") +
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
    panel.background = element_rect(fill = "#a9bf99"),
    legend.position = "top",
    panel.border = element_rect(fill = NA))

## Save plot ----
ggsave(
  "output/plots/map_potts.png",
  width = 0.48 * 7.5,
  height = 6,
  units = "in",
  dpi = "retina")



# Janky ecozone map (ignore) ----------------------------------------------

zonemap <- function(plot) {
  ggplot() + 
    geom_spatraster(
      data = 
        ecozone_rast |> 
        crop(
          plots |> 
          filter(plot_name == plot) |> 
          st_transform(
            crs = st_crs(ecozone_rast)) |> 
            vect()),
      maxcell = Inf) + 
    scale_fill_manual(
    values = c(
      "gray80", "forestgreen", "violet", "lightgreen", 
      "lightgoldenrod1", "#0978AB", "turquoise3", "white"),
    labels = c(
      "barrens", "cove", "dry oak-heath", "dry pine-oak", 
      "dry-mesic oak", "floodplains", "mesic oak", "unknown"), 
    drop = FALSE) +
    labs(
      title = paste0("\n", plot),
      subtitle = 
        str_replace(filter(plots, plot_name == plot)$plot_description, "Road", "Rd."),
      fill = "Ecozone") +
    theme_void() +
    theme(
      legend.position = "bottom",
      plot.subtitle = element_text(size = 8),
      panel.border = element_rect(
        fill = NA, color = "black", linewidth = 0.75)) +
    geom_sf(
      data = points |> filter(plot_name == plot), 
      shape = 1, size = 3, stroke = 1.25) +
    coord_sf(expand = F)
}

plot_list <- list()
plots_to_plot <- 
  plots |> filter(plot_type == "Population") |> pull(plot_name) |> unique()

for(i in 1:length(plots_to_plot)) {
  plot_list[[i]] <- zonemap(plots_to_plot[[i]])}

library(ggpubr)

arranged_plots <- 
  ggarrange(
    plot_list[[1]], 
    plot_list[[2]],
    plot_list[[3]],
    plot_list[[4]],
    plot_list[[5]],
    plot_list[[6]],
    plot_list[[7]],
    plot_list[[8]],
    plot_list[[9]],
    plot_list[[10]],
    plot_list[[11]],
    plot_list[[12]],
    plot_list[[13]],
    plot_list[[14]],
    plot_list[[15]],
    nrow = 5, ncol = 3,
    common.legend = T)

arranged_plots
ggsave(
  filename = "output/plots/ecozone_comboplots.png",
  width = 6, height = 10, units = "in", dpi = "retina")


