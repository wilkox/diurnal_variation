# Libraries
library(tidyverse)
library(fs)
library(furrr)

# Set up environment
Sys.setenv(SOURCETRACKER_PATH = "../opt/sourcetracker-1.0.1/")
plan("multiprocess")

# Create output directory
if (dir_exists("source_predictions")) dir_delete("source_predictions")
dir_create("source_predictions")

# Prepare a list of runs
runs <- dir_ls("maps") %>%
  enframe(value = "map_path") %>%
  select(-name) %>%
  mutate(map_path = path_wd(map_path)) %>%
  mutate(run = map_path %>% path_file() %>% path_ext_remove() %>% as.character()) %>%
  mutate(OTU_table_path = path_wd("OTU_tables", run, ext = "tsv")) %>%
  filter(file_exists(OTU_table_path)) %>%
  mutate(output_dir = path_wd("source_predictions", run)) %>%
  select(run, everything())

# Run SourceTracker
sourcetracker_path <- "../opt/sourcetracker-1.0.1/sourcetracker_for_qiime.r"
run_sourcetracker <- function(OTU_table_path, map_path, output_dir, ...) {
  message("Running for ", OTU_table_path)
  system2(
    "Rscript",
    c(sourcetracker_path, "-i", OTU_table_path, "-m", map_path, "-o", output_dir)
  )
}

future_pmap(runs, run_sourcetracker)
