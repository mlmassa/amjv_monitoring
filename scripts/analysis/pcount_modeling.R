
# Setup -------------------------------------------------------------------

# Load required packages
library(tidyverse)
library(unmarked)
library(MuMIn)

# Allow tidy parsing of unmarked results
source("scripts/wrangling/tidy_method_unmarked.R")

# Load unmarked frames (created in format_umf.R)
read_rds("data/processed/umfs.rds") |> 
  list2env(.GlobalEnv)

# Set focal species

# umf <- WOTH
# sp_long <- "Wood Thrush"
# sp <- "WOTH"

# umf <- SCTA
# sp_long <- "Scarlet Tanager"
# sp <- "SCTA"

umf <- WEWA
sp_long <- "Worm-eating Warbler"
sp <- "WEWA"

# List available variables
# names(umf@obsCovs)
# names(umf@siteCovs)

# What are we missing?
# tibble(umf@siteCovs) |> 
#   map(~sum(is.na(.x))) |> 
#   bind_cols()

# Can't use basal-- too many missing values
  
# Check collinearity ------------------------------------------------------
# Observation covs that should not go in a model together
# umf@obsCovs |> 
#   select(where(is.numeric)) |> 
#   cor(use = "pairwise.complete.obs") |> 
#   as.data.frame() |> 
#   rownames_to_column("var1") |> 
#   pivot_longer(-1, names_to = "var2", values_to = "cor") |> 
#   filter(cor != 1, abs(cor) > 0.5) |> 
#   group_by(cor) |> slice_head(n = 1)

# Site covariates that should not go in a model together
# umf@siteCovs |> 
#   select(where(is.numeric)) |> 
#   cor(use = "pairwise.complete.obs") |> 
#   as.data.frame() |> 
#   rownames_to_column("var1") |> 
#   pivot_longer(-1, names_to = "var2", values_to = "cor") |> 
#   filter(cor != 1, abs(cor) > 0.5) |> 
#   group_by(cor) |> slice_head(n = 1)

# Basal area and canopy height are slightly related anyway, so let's ditch basal

# We are using K=40 because there is unlikely to be an abundance higher than 50, although it is mathematically possible in a model. It also cuts down on calculation time. For discussion of this, see https://groups.google.com/g/unmarked/c/KaU4nyxhe5E

# Set up parallel ---------------------------------------------------------
library(parallel)

clust <- 
  makeCluster(getOption("cl.cores", 4), type = "PSOCK")

clusterEvalQ(clust, library(unmarked))
clusterExport(clust, "umf")

# stats4 is needed for AIC to work with unmarkedFit objects
library(stats4)
clusterCall(clust, "library", "stats4", character.only = TRUE)

# Detection ---------------------------------------------------------------

# No temperature (single year)
# No observer (fails to converge and no difference)

m.d.global <-
  pcount(
    data = umf,
    formula = 
      ~ wind + scale(cloud) + precip + scale(start_sun) + scale(doy)
      ~ 1,
    K = 40,
    mixture = "ZIP");beepr::beep("treasure")

system.time(
  m.d.dredge <-
    dredge(
      m.d.global,
      cluster = clust,
      rank = AIC,
      ));beepr::beep("treasure")

d.formula <-
  get.models(m.d.dredge, subset = 1)[[1]]@formula[[2]] |> 
  as.character() |> 
  pluck(2)

#filter(m.d.dredge, delta < 2)

#message(sp, " detection: ", d.formula)

# *SCTA: cloud + start_sun + wind
# *WEWA: cloud + start_sun + doy + wind
# *WOTH: cloud + start_sun

# Abundance ---------------------------------------------------------------

# No basal (missing data)
# 15 seconds
m.a.global <-
  pcount(
    data = umf,
    formula = 
      paste(
        "~", d.formula,
        "~ scale(for_1500) + scale(ed_buff) + scale(for_5000) + scale(d_stream) + scale(canopy_h) + scale(canopy_cv) + year") |> 
      as.formula(),
    K = 40,
    mixture = "ZIP");beepr::beep("treasure")

# Ecozone only
m.a.eco <-
  pcount(
    data = umf,
    formula = 
      paste(
        "~", d.formula,
        "~ ecozone") |> 
      as.formula(),
    K = 40,
    mixture = "ZIP");beepr::beep("treasure")

# Ownership only
m.a.own <-
  pcount(
    data = umf,
    formula = 
      paste(
        "~", d.formula,
        "~ ownership") |> 
      as.formula(),
    K = 40,
    mixture = "ZIP");beepr::beep("treasure")

# Basal only
m.a.basal <-
  pcount(
    data = umf,
    formula = 
      paste(
        "~", d.formula,
        "~ scale(basal)") |> 
      as.formula(),
    K = 40,
    mixture = "ZIP");beepr::beep("treasure")

# Clear memory
gc()

# Set model subsets:
# Set detection model terms (MUST be included, differ by species, see above)
# Disallow correlated terms (Same for all species)
if(sp == "WOTH") {
  subset <- expression(
    `p(scale(cloud))` && `p(scale(start_sun))` 
    && !(`lam(scale(for_1500))` & `lam(scale(ed_buff))`) 
    && !(`lam(scale(canopy_h))` & `lam(scale(canopy_cv))`))
} else if (sp == "SCTA") {
  subset <- expression(
    `p(scale(cloud))` && `p(scale(start_sun))` && `p(wind)` 
    && !(`lam(scale(for_1500))` & `lam(scale(ed_buff))`) 
    && !(`lam(scale(canopy_h))` & `lam(scale(canopy_cv))`))
} else if (sp == "WEWA") {
  subset <- expression(
    `p(scale(cloud))` && `p(scale(start_sun))` && `p(wind)` && `p(scale(doy))` 
    && !(`lam(scale(for_1500))` & `lam(scale(ed_buff))`) 
    && !(`lam(scale(canopy_h))` & `lam(scale(canopy_cv))`))
} else print("Try again")

system.time(
  m.a.dredge <-
    dredge(
      m.a.global,
      cluster = clust,
      rank = AIC,
      subset = subset
      ));beepr::beep("treasure")

# Best models
#filter(m.a.dredge, delta < 2)

write_rds(m.a.dredge, paste0("output/aictab_", sp, ".rds"))
# Top model
# Get best model
m.best <-
  get.models(m.a.dredge, subset = 1)[[1]]

write_rds(m.best, paste0("output/best_model_", sp, ".rds"))

a.formula <-
  m.best@formula


# Close parallel ----------------------------------------------------------

stopCluster(clust)

# Generate predictions ----------------------------------------------------

m.best <- read_rds(paste0("output/best_model_", sp, ".rds"))

# Set defaults

pred_defaults <-
  data.frame(
    wind = umf@obsCovs$wind[1],
    precip = umf@obsCovs$precip[5],
    cloud = mean(umf@obsCovs$cloud, na.rm = T),
    doy = mean(umf@obsCovs$doy, na.rm = T),
    start_sun = mean(umf@obsCovs$start_sun, na.rm = T))

pred_abund_defaults <-
  data.frame(
    basal = mean(umf@siteCovs$basal, na.rm = T),
    for_1500 = mean(umf@siteCovs$for_1500, na.rm = T),
    for_5000 = mean(umf@siteCovs$for_5000, na.rm = T),
    year = umf@siteCovs$year[1], 
    para = mean(umf@siteCovs$para, na.rm = T),
    canopy_h = mean(umf@siteCovs$canopy_h, na.rm = T),
    canopy_cv = mean(umf@siteCovs$canopy_cv, na.rm = T),
    core = mean(umf@siteCovs$core, na.rm = T),
    d_stream = mean(umf@siteCovs$d_stream, na.rm = T),
    ed_buff = mean(umf@siteCovs$ed_buff, na.rm = T))



# Create predictions
preds <- list(
  # Detection covariates
  wind = 
    predict(m.best,
    newdata = bind_cols(
      wind = levels(umf@obsCovs$wind),
      select(pred_defaults, -wind)),
    type = "det", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  precip = 
    predict(m.best,
    newdata = bind_cols(
      precip = levels(umf@obsCovs$precip),
      select(pred_defaults, -precip)),
    type = "det", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  cloud = 
    predict(m.best,
    newdata = bind_cols(
      cloud = seq(
        min(umf@obsCovs$cloud, na.rm = T), 
        max(umf@obsCovs$cloud, na.rm = T), 
        length.out = 20),
      select(pred_defaults, - cloud)),
    type = "det", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  doy = 
    predict(m.best, 
    newdata = bind_cols(
      select(pred_defaults, -doy),
      doy = seq(
        min(umf@obsCovs$doy, na.rm = T), 
        max(umf@obsCovs$doy, na.rm = T), 
        length.out = 20)),
    type = "det", appendData = T) |> 
    mutate(
      doy = as.Date(doy - 1, origin = "2023-01-01"),
      sp = sp, 
      sp_long = sp_long),
  start_sun = 
    predict(m.best,
    newdata = bind_cols(
      start_sun = seq(
        min(umf@obsCovs$start_sun, na.rm = T), 
        max(umf@obsCovs$start_sun, na.rm = T), 
        length.out = 20),
      select(pred_defaults, - start_sun)),
    type = "det", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  # Solo abundance models
  ownership = 
    predict(m.a.own,
    newdata = bind_cols(ownership = levels(umf@siteCovs$ownership)),
    type = "state", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  ecozone = 
    predict(m.a.eco,
    newdata = bind_cols(ecozone = levels(umf@siteCovs$ecozone)),
    type = "state", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  basal =
    predict(m.a.basal,
    newdata = bind_cols(
      basal = seq(
        min(umf@siteCovs$basal, na.rm = T), 
        max(umf@siteCovs$basal, na.rm = T), 
        length.out = 20)),
    type = "state", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  # Abundance covariates
  year = 
    predict(m.best,
    newdata = 
      bind_cols(
        year = levels(umf@siteCovs$year),
        select(pred_abund_defaults, -year)),
    type = "state", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  for_1500 =
    predict(m.best,
    newdata = 
      bind_cols(
        for_1500 = seq(
          min(umf@siteCovs$for_1500, na.rm = T), 
          max(umf@siteCovs$for_1500, na.rm = T), 
          length.out = 20),
        select(pred_abund_defaults, -for_1500)),
    type = "state", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  for_5000 = 
    predict(m.best,
    newdata = 
      bind_cols(
        for_5000 = seq(
          min(umf@siteCovs$for_5000, na.rm = T), 
          max(umf@siteCovs$for_5000, na.rm = T), 
          length.out = 20),
        select(pred_abund_defaults, -for_5000)),
    type = "state", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  core = 
    predict(m.best,
    newdata = 
      bind_cols(
        core = seq(
          min(umf@siteCovs$core, na.rm = T), 
          max(umf@siteCovs$core, na.rm = T), 
          length.out = 20),
        select(pred_abund_defaults, -core)),
    type = "state", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  para = 
    predict(m.best,
    newdata = 
      bind_cols(
        para = seq(
          min(umf@siteCovs$para, na.rm = T), 
          max(umf@siteCovs$para, na.rm = T), 
          length.out = 20),
        select(pred_abund_defaults, -para)),
    type = "state", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  d_stream = 
    predict(m.best,
    newdata = 
      bind_cols(
        d_stream = seq(
          min(umf@siteCovs$d_stream, na.rm = T), 
          max(umf@siteCovs$d_stream, na.rm = T), 
          length.out = 20),
        select(pred_abund_defaults, -d_stream)),
    type = "state", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  ed_buff = 
    predict(m.best,
    newdata = 
      bind_cols(
        ed_buff = seq(
          min(umf@siteCovs$ed_buff, na.rm = T), 
          max(umf@siteCovs$ed_buff, na.rm = T), 
          length.out = 20),
        select(pred_abund_defaults, -ed_buff)),
    type = "state", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  canopy_h = 
    predict(m.best,
    newdata = 
      bind_cols(
        canopy_h = seq(
          min(umf@siteCovs$canopy_h, na.rm = T), 
          max(umf@siteCovs$canopy_h, na.rm = T), 
          length.out = 20),
        select(pred_abund_defaults, -canopy_h)),
    type = "state", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long),
  canopy_cv = 
    predict(m.best,
    newdata = 
      bind_cols(
        canopy_cv = seq(
          min(umf@siteCovs$canopy_cv, na.rm = T), 
          max(umf@siteCovs$canopy_cv, na.rm = T), 
          length.out = 20),
        select(pred_abund_defaults, -canopy_cv)),
    type = "state", appendData = T) |> 
    mutate(sp = sp, sp_long = sp_long)
  )

# Write data
write_rds(preds, paste0("output/preds_", sp, ".rds"))

# Generate model-averaged predictions -------------------------------------

system.time(
  m.avg <-
    model.avg(
      m.a.dredge, 
      subset = delta < 2, 
      fit = T));beepr::beep("treasure")

# Write model
write_rds(m.avg, paste0("output/model_avg_", sp, ".rds"))

m.avg <- read_rds(paste0("output/model_avg_", sp, ".rds"))

## Pred avg list ----
preds_avg <- list(
  # Detection covariates
  wind = 
    predict(
      m.avg,
      newdata = bind_cols(
        wind = levels(umf@obsCovs$wind),
        select(pred_defaults, -wind)),
      type = "det") |> 
    bind_cols() |> 
    mutate(
      wind = levels(umf@obsCovs$wind),
      sp = sp, sp_long = sp_long),
  precip = 
    predict(
      m.avg,
      newdata = bind_cols(
        precip = levels(umf@obsCovs$precip),
        select(pred_defaults, -precip)),
      type = "det") |> 
    bind_cols() |> 
    mutate(
      precip = levels(umf@obsCovs$precip),
      sp = sp, sp_long = sp_long),
  cloud = 
    predict(
      m.avg,
      newdata = bind_cols(
        cloud = seq(
          min(umf@obsCovs$cloud, na.rm = T), 
          max(umf@obsCovs$cloud, na.rm = T), 
          length.out = 20),
        select(pred_defaults, - cloud)),
      type = "det") |> 
    bind_cols() |> 
    mutate(
      cloud = seq(
        min(umf@obsCovs$cloud, na.rm = T), 
        max(umf@obsCovs$cloud, na.rm = T), 
        length.out = 20),
      sp = sp, sp_long = sp_long),
  doy = 
    predict(
      m.avg, 
      newdata = bind_cols(
        select(pred_defaults, -doy),
        doy = seq(
          min(umf@obsCovs$doy, na.rm = T), 
          max(umf@obsCovs$doy, na.rm = T), 
          length.out = 20)),
      type = "det") |> 
    bind_cols() |> 
    mutate(
      doy = seq(
        min(umf@obsCovs$doy, na.rm = T), 
        max(umf@obsCovs$doy, na.rm = T), 
        length.out = 20),
      doy = as.Date(doy - 1, origin = "2023-01-01"),
      sp = sp, sp_long = sp_long),
  start_sun = 
    predict(
      m.avg,
      newdata = bind_cols(
        start_sun = seq(
          min(umf@obsCovs$start_sun, na.rm = T), 
          max(umf@obsCovs$start_sun, na.rm = T), 
          length.out = 20),
        select(pred_defaults, - start_sun)),
      type = "det") |> 
    bind_cols() |> 
    mutate(
      start_sun = seq(
        min(umf@obsCovs$start_sun, na.rm = T), 
        max(umf@obsCovs$start_sun, na.rm = T), 
        length.out = 20),
      sp = sp, sp_long = sp_long),
  # Abundance covariates
  year = 
    predict(
      m.avg,
      newdata = 
        bind_cols(
          year = levels(umf@siteCovs$year),
          select(pred_abund_defaults, -year)),
      type = "state") |> 
    bind_cols() |> 
    mutate(
      year = levels(umf@siteCovs$year),
      sp = sp, sp_long = sp_long),
  for_1500 =
    predict(
      m.avg,
      newdata = 
        bind_cols(
          for_1500 = seq(
            min(umf@siteCovs$for_1500, na.rm = T), 
            max(umf@siteCovs$for_1500, na.rm = T), 
            length.out = 20),
          select(pred_abund_defaults, -for_1500)),
      type = "state") |> 
    bind_cols() |> 
    mutate(
      for_1500 = seq(
            min(umf@siteCovs$for_1500, na.rm = T), 
            max(umf@siteCovs$for_1500, na.rm = T), 
            length.out = 20),
      sp = sp, sp_long = sp_long),
  for_5000 = 
    predict(
      m.avg,
      newdata = 
        bind_cols(
          for_5000 = seq(
            min(umf@siteCovs$for_5000, na.rm = T), 
            max(umf@siteCovs$for_5000, na.rm = T), 
            length.out = 20),
          select(pred_abund_defaults, -for_5000)),
      type = "state") |> 
    bind_cols() |> 
    mutate(
      for_5000 = seq(
            min(umf@siteCovs$for_5000, na.rm = T), 
            max(umf@siteCovs$for_5000, na.rm = T), 
            length.out = 20),
      sp = sp, sp_long = sp_long),
  core = 
    predict(
      m.avg,
      newdata = 
        bind_cols(
          core = seq(
            min(umf@siteCovs$core, na.rm = T), 
            max(umf@siteCovs$core, na.rm = T), 
            length.out = 20),
          select(pred_abund_defaults, -core)),
    type = "state") |> 
    bind_cols() |> 
    mutate(
      core = seq(
            min(umf@siteCovs$core, na.rm = T), 
            max(umf@siteCovs$core, na.rm = T), 
            length.out = 20),
      sp = sp, sp_long = sp_long),
  para = 
    predict(
      m.avg,
      newdata = 
        bind_cols(
          para = seq(
            min(umf@siteCovs$para, na.rm = T), 
            max(umf@siteCovs$para, na.rm = T), 
            length.out = 20),
          select(pred_abund_defaults, -para)),
      type = "state") |> 
    bind_cols() |> 
    mutate(
      para = seq(
            min(umf@siteCovs$para, na.rm = T), 
            max(umf@siteCovs$para, na.rm = T), 
            length.out = 20),
      sp = sp, sp_long = sp_long),
  d_stream = 
    predict(
      m.avg,
      newdata = 
        bind_cols(
          d_stream = seq(
            min(umf@siteCovs$d_stream, na.rm = T), 
            max(umf@siteCovs$d_stream, na.rm = T), 
            length.out = 20),
          select(pred_abund_defaults, -d_stream)),
      type = "state") |> 
    bind_cols() |> 
    mutate(
      d_stream = seq(
            min(umf@siteCovs$d_stream, na.rm = T), 
            max(umf@siteCovs$d_stream, na.rm = T), 
            length.out = 20),
      sp = sp, sp_long = sp_long),
  ed_buff = 
    predict(
      m.avg,
      newdata = 
        bind_cols(
          ed_buff = seq(
            min(umf@siteCovs$ed_buff, na.rm = T), 
            max(umf@siteCovs$ed_buff, na.rm = T), 
            length.out = 20),
          select(pred_abund_defaults, -ed_buff)),
      type = "state") |> 
    bind_cols() |> 
    mutate(
      ed_buff = seq(
            min(umf@siteCovs$ed_buff, na.rm = T), 
            max(umf@siteCovs$ed_buff, na.rm = T), 
            length.out = 20),
      sp = sp, sp_long = sp_long),
  canopy_h = 
    predict(
      m.avg,
      newdata = 
        bind_cols(
          canopy_h = seq(
            min(umf@siteCovs$canopy_h, na.rm = T), 
            max(umf@siteCovs$canopy_h, na.rm = T), 
            length.out = 20),
          select(pred_abund_defaults, -canopy_h)),
      type = "state") |> 
    bind_cols() |> 
    mutate(
      canopy_h = seq(
            min(umf@siteCovs$canopy_h, na.rm = T), 
            max(umf@siteCovs$canopy_h, na.rm = T), 
            length.out = 20),
      sp = sp, sp_long = sp_long),
  canopy_cv = 
    predict(
      m.avg,
      newdata = 
        bind_cols(
          canopy_cv = seq(
            min(umf@siteCovs$canopy_cv, na.rm = T), 
            max(umf@siteCovs$canopy_cv, na.rm = T), 
            length.out = 20),
          select(pred_abund_defaults, -canopy_cv)),
      type = "state") |> 
    bind_cols() |> 
    mutate(
      canopy_cv = seq(
            min(umf@siteCovs$canopy_cv, na.rm = T), 
            max(umf@siteCovs$canopy_cv, na.rm = T), 
            length.out = 20),
      sp = sp, sp_long = sp_long)
  )

# Write data
write_rds(preds_avg, paste0("output/preds_avg_", sp, ".rds"))
