# Libraries
library(tidyverse)
library(furrr)
plan(multisession)

# Function to decompose a time series of species abundances
decompose_abundances <- function(x) {

  x %>%
    select(timepoint, abundance) %>%
    mutate_at("timepoint", ~ factor(., levels = min(.):max(.))) %>%
    complete(timepoint, fill = list(abundance = NA)) %>%
    arrange(timepoint) %>%
    pull(abundance) %>%
    ts(frequency = 2)  %>%
    imputeTS::na_seadec() %>%
    decompose("additive") %>%
    .[1:4] %>%
    map(as.double) %>%
    as_tibble() %>%
    mutate(timepoint = 1:n()) %>%
    select(timepoint, observed = x, everything())
}

# Function to calculate seasonality strength
seasonality_strength <- function(time_series) {
  time_series <- slice(time_series, -1, -nrow(time_series))
  ss <- 1 - (sd(time_series$random) / sd(time_series$random + time_series$seasonal))
  ss <- replace_na(ss, 0)
  max(c(0, ss))
}

# Function to calculate permuted seasonality strength
permuted_seasonality_strengths <- function(x, permutations = 999) {
  # Generate random permutations of abundance
  tibble(
    permutation = 1:permutations,
    abundances = rep(list(select(x, timepoint, abundance)), permutations)
    ) %>%
    mutate(abundances = map(abundances, 
                            ~ mutate(.x, abundance = sample(abundance)))) %>%
  mutate(ts = map(abundances, decompose_abundances)) %>%
  mutate(ss = map_dbl(ts, seasonality_strength)) %>%
  pull(ss)
}

# Load list of species that appear in negative control samples
nc_species <- "../abundances/control_species_abundances.tsv" %>%
  read_tsv() %>%
  filter(abundance > 0) %>%
  pull(species) %>%
  unique()

# Load abundances, remove the negative control species, and renormalise
# remaining species abundances
abundances <- read_tsv("../abundances/species_abundances_decontaminated.tsv") %>%
  select(sample, site, location, timepoint, species, abundance) %>%
  filter(! species %in% nc_species) %>%
  group_by(sample, site, location, timepoint) %>%
  mutate(abundance = 100 * abundance / sum(abundance)) %>%
  ungroup()

# Run time series decompositions and seasonality strength calculations
seasonality_strengths <- abundances %>%
  select(site, location, timepoint, species, abundance) %>%
  add_count(site, location, species, wt = abundance, name = "total_abundance") %>%
  filter(total_abundance > 0) %>%
  select(-total_abundance) %>% 
  nest(abundances = c(timepoint, abundance)) %>%
  mutate(ts = map(abundances, decompose_abundances)) %>%
  mutate(observed_ss = map_dbl(ts, seasonality_strength)) %>%
  select(-ts) %>%
  mutate(permuted_ss = future_map(abundances, permuted_seasonality_strengths,
         .progress = TRUE)) %>%
  select(-abundances)

# Calculate p values
seasonality_strengths <- seasonality_strengths %>%
  mutate(p = map2_dbl(observed_ss, permuted_ss, ~ sum(.y >= .x) / length(.y))) %>%
  select(-permuted_ss) %>%
  rename(seasonality_strength = observed_ss)

# Write to file
write_tsv(seasonality_strengths, "no_nc_seasonality_strengths.tsv")
