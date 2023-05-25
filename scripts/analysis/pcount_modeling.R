
# Setup -------------------------------------------------------------------

# Load required packages
library(tidyverse)
library(unmarked)
library(AICcmodavg)
library(MuMIn)

# Allow tidy parsing of unmarked results
source("scripts/wrangling/tidy_method_unmarked.R")

# Load unmarked frames (created in format_umf.R)
read_rds("data/processed/umfs.rds") |> 
  list2env(.GlobalEnv)

umf <- SCTA
sp_long <- "Scarlet Tanager"
sp <- "SCTA"

# Set theme
theme_set(
  theme_bw())

# What variables are available?
str(umf)

# Detection ---------------------------------------------------------------

## One by one ----


m.d1 <-
  pcount(
    data = umf,
    formula = ~ start_sun ~ 1,
    K = 50,
    mixture = "ZIP")

m.d2 <-
  pcount(
    data = umf,
    formula = ~ doy ~ 1,
    K = 50,
    mixture = "ZIP")

predict(
  m.d1,
  newdata = data.frame(
    start_sun = 
      seq(
        min(umf@obsCovs$start_sun, na.rm = T), 
        max(umf@obsCovs$start_sun, na.rm = T), 
        by = 10)),
  type = "det",
  appendData = T) |> 
  ggplot(aes(x = start_sun))+
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.3) +
  geom_line(aes(y = Predicted)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 60, 120, 180, 240, 300)) +
  labs(
    y = "Predicted detection",
    x = "Minutes since sunrise",
    title = sp_long) +
  coord_cartesian(expand = T, clip = "off")

predict(
  m.d2,
  newdata = data.frame(
    doy = 
      seq(
        min(umf@obsCovs$doy, na.rm = T), 
        max(umf@obsCovs$doy, na.rm = T), 
        by = 10)),
  type = "det",
  appendData = T) |> 
  mutate(doy = as.Date(doy - 1, origin = "2023-01-01")) |> 
  ggplot(aes(x = doy))+
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.3) +
  geom_line(aes(y = Predicted)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_date(
    date_breaks = "7 days",
    labels = scales::label_date("%b %d")) +
  labs(
    y = "Predicted detection",
    x = "Day of year",
    title = sp_long) +
  coord_cartesian(expand = T, clip = "off")


## Dredge ----

m.d.global <-
  pcount(
    data = umf,
    formula = ~ wind + cloud + precip + start_sun + doy ~ 1,
   # K = 50,
    mixture = "ZIP")

dredge(
  m.d.global)

# WOTH: start_sun + cloud + wind

pcount(
    data = umf,
    formula = ~ scale(start_sun) ~ scale(basal) + year,
   # K = 50,
    mixture = "ZIP")
  