# This script imports data. Run it prior to any major analysis.

# Setup -------------------------------------------------------------------

# Load required packages
library(googledrive)
library(googlesheets4)
library(tidyverse)

# Authorize Google Drive access
# You may need to obtain a new token then run again.
drive_auth()

# Authorize Google Sheets access
gs4_auth(token = drive_token())

# Import data from Google Sheets ------------------------------------------

# Get vector of sheet names
sheets <-
  drive_get('AMJV_data_entry_tidy') %>% 
  gs4_get() %>% 
  sheet_names()

# Read in data. This will take a minute.
sheets %>% 
    map(
    function(name) {
      drive_get('AMJV_data_entry_tidy') %>% 
      gs4_get() %>% 
      read_sheet(name)
    }) %>% 
    set_names(sheets) %>% 
# Import tables separately to global environment
list2env(.GlobalEnv)

# Assign eventID ----------------------------------------------------------

# Function: random string generator ABC1234. 175,760,000 options.
random_string <-
  function(n = 5000) {
    a <- 
      do.call(paste0, replicate(3, sample(LETTERS, n, TRUE), FALSE))
    paste0(
      a, 
      sprintf('%04d', sample(9999, n, TRUE)))
  }

# Generate a random string for all event_ids
id_strings <-
  random_string(n = length(count_data$event_id)) 

# Assign id strings to event_ids
count_data$event_id <- id_strings

# Clean up
rm(id_strings, random_string)

# Write data --------------------------------------------------------------

# Create list of tibbles
amjv_data <-
  mget(sheets)

# Save objects as RDS.
# This data will be loaded for other operations.
# Rerun THIS SCRIPT whenever you do analysis if the database has been edited.
write_rds(
  amjv_data,
  'data/raw/amjv_data.rds')
  