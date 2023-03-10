# This script is intended to generate a list of all species that occur in
# monitored counties during the survey period. It can be rerun with different
# locations or different date windows. Use the resultant species list as
# a taxonomy table and for data validation.

#  Setup ------------------------------------------------------------------

library(tidyverse)
library(rebird)

# Get eBird frequency data ------------------------------------------------

# Bootleg version of rebird::ebirdfreq()
# Function is broken because you need to log in to eBird to access frequency
bootleg_ebirdfreq <-
  function (
    loctype, 
    loc, 
    startyear = 1900, 
    endyear = format(Sys.Date(), "%Y"), 
    startmonth = 1, 
    endmonth = 12, 
    long = TRUE, ...) {
    # Ensure month inputs are valid
    if (!all(c(startmonth, endmonth) %in% 1:12)) {
      stop("Invalid month(s) provided. Must be integer between 1 and 12.")
    }
    # Get current year
    currentyear <- format(Sys.Date(), "%Y")
    # Ensure year inputs are valid
    if (!all(c(startyear, endyear) %in% 1900:currentyear)) {
      stop(paste0("Invalid year(s) provided. Must be integer between 1900 and ", 
        currentyear, "."))
    }
    # Ensure hotspot input is valid
    if (loctype == "hotspots") {
      if (!grepl("^L\\d{1,8}$", loc)) 
        stop("Invalid hotspot code")
    }
    else if (loctype %in% c("counties", "states")) {
    }
    else {
      stop("Not a valid location type")
    }
    args <- list(r = loc, bmo = startmonth, emo = endmonth, 
      byr = startyear, eyr = endyear, personal = "false", 
      fmt = "tsv")
    url <- "http://ebird.org/ebird/barchartData"
    url_full <- httr::modify_url(url, query = args)
    message(url_full)
  }

# List locations to base expected frequency on
my_locations <-
  c(
    # Bath, VA (surveyed):
    "US-VA-017",
    # Rockbridge, VA (surveyed):
    "US-VA-163",
    # Highland, VA (nearby):
    "US-VA-091",
    # Augusta, VA (nearby):
    "US-VA-015",
    # Alleghany, VA (nearby):
    "US-VA-005",
    # Pocahontas, WV (nearby):
    "US-WV-075",
    # Pendleton, WV (nearby):
    "US-WV-071")

# Run bootleg URL generator over all locations:
mapply(
  bootleg_ebirdfreq,
  loctype = "counties",
  loc = my_locations,
  startmonth = 4,
  endmonth = 7)

# !!!IMPORTANT!!!

# Open the message URLs in browser to download files.
# Move the files to data/raw before proceeding.
# Automating this process is not possible with current eBird API rules.

# Get frequencies ---------------------------------------------------------

# Create container. We'll start with WOTH because we know they're present!
species_list <- "Wood Thrush"

# For loop: generate species lists and append them to species_list
for(i in my_locations) {
  
  # Read in frequency table:
    freq <- 
    read.delim(
      paste0(
        "data/raw/ebird_",
        i,
        "__1900_2021_4_7_barchart.txt"), 
      skip = 12, 
      stringsAsFactors = FALSE) %>% 
    as_tibble() %>% 
    select(-50) %>% 
    # Rename columns: name and "Month-Week#"
    magrittr::set_colnames(
      c('name', 
        vapply(
          month.abb, 
          paste, 
          FUN.VALUE = character(4), 
          1:4, 
          sep = "-"))) 
  
  # Transform into long data with clean names:
    freq_long <-
    # Get long data
    reshape(
      data.frame(freq[-1, ]), 
      varying = 2:49, 
      direction = "long", 
      v.names = "frequency", 
      idvar = "name", 
      timevar = "month_week", 
      times = names(freq)[2:49]) %>% 
      tbl_df() %>% 
    # Remove HTML if present (if you got scientific names)
    mutate_at("name", str_replace, "\\s\\(<em\\sclass=sci>", "") %>% 
    mutate_at("name", str_replace, "<\\Wem>\\)", "") %>% 
    # Filter the data
    filter(
      # Only May and June
      month_week %in% c("May-1", "May-2", "May-3", "May-4", 
                        "Jun-1", "Jun-2", "Jun-3", "Jun-4"),
      # Species present
      frequency != 0)
  
  # Create list of species for that location 
    new_species <- 
    unique(freq_long$name)
  
  # Append to full species list, no duplicates
    for(s in new_species) {
  species_list[[length(species_list) + 1]] <- s
  species_list <- unique(species_list)
  }
  
  # Clean up
    rm(freq, freq_long, new_species)
  
}

# Get taxonomy ------------------------------------------------------------

# Update taxonomy (DCCO recently changed sci_name. This may take a minute.)
new_tax <- ebirdtaxonomy()

# Generate species table
species <-
  tibble(comName = species_list) %>% 
    # Add eBird taxonomy
    left_join(
      new_tax,
      by = "comName") %>% 
    select(
      sciName,
      comName,
      category,
      taxonOrder,
      bandingCodes) %>% 
  filter(
    !is.na(bandingCodes)) %>% 
  arrange(taxonOrder)

# Write csv of species
write_csv(
  species,
  "data/processed/species_list.csv")
