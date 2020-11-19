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

# Load abundances and calculate seasonality strengths
# I'll do this in three separate tranches to prevent the process being killed
# for memory pressure
do_tranche <- function(range) {
  seasonality_strengths <- read_tsv("../abundances/species_abundances_decontaminated.tsv") %>%
    select(site, location, timepoint, species, abundance) %>%
    add_count(site, location, species, wt = abundance, name = "total_abundance") %>%
    filter(total_abundance > 0) %>%
    select(-total_abundance) %>% 
    nest(abundances = c(timepoint, abundance)) %>%
    slice(range) %>%
    mutate(ts = map(abundances, decompose_abundances)) %>%
    mutate(observed_ss = map_dbl(ts, seasonality_strength)) %>%
    select(-ts) %>%
    mutate(permuted_ss = future_map(abundances, permuted_seasonality_strengths,
           .progress = TRUE)) %>%
    select(-abundances)
}

t1 <- do_tranche(0:800)
t2 <- do_tranche(801:1600)
t3 <- do_tranche(1601:2400)
seasonality_strengths <- bind_rows(t1, t2, t3)

# Calculate p values
seasonality_strengths <- seasonality_strengths %>%
  mutate(p = map2_dbl(observed_ss, permuted_ss, ~ sum(.y >= .x) / length(.y))) %>%
  select(-permuted_ss) %>%
  rename(seasonality_strength = observed_ss)

# Write to file
write_tsv(seasonality_strengths, "seasonality_strengths.tsv")
