library(tidyverse)
library(WorldFlora)
library(doParallel)


# Load sPlot and xylem data
load("data/sPlot/Xylem_sPlot.RData")
traits <- readRDS("data/cleaned/trait_data.rds") %>%
  rename(species = species_wfo)


# merge trait and sPlot data

d <- DT.xylem %>% 
  rename_all(tolower) %>%
  distinct(plot_id = plotobservationid, species)

plots <- header.xylem %>%
  janitor::clean_names() %>%
  rename(plot_id = 1) %>%
  select(plot_id, longitude, latitude, location_uncertainty_m, continent, s_biome, ecoregion) %>%
  filter(location_uncertainty_m <= 500) %>%
  rename(biome = s_biome) %>%
  na.omit() %>%
  distinct()


# merge with plot data and drop ecoregions with fewer than 100 observations
d2 <- inner_join(traits, d) %>% 
  inner_join(plots) %>%
  add_count(ecoregion) %>%
  filter(n > 100) %>%
  select(plot_id, longitude, latitude, location_uncertainty_m, continent, biome, ecoregion, species, p50, rd_mean, rd_max)


# prepare merged sPlot-trait data
out <- d2 %>%
  select(plot_id, species, p50, rd_max) %>%
  arrange(plot_id, species)

# prepare sPlot presence data
pres <- out %>%
  select(plot_id, species) %>%
  distinct()

# prepare trait data
trait2 <- d2 %>%
  select(species, p50, rd_max) %>%
  distinct() %>%
  arrange(species)

# prepare plot data
plots2 <- d2 %>%
  select(plot_id:ecoregion) %>%
  distinct()

# export data
saveRDS(plots2, "data/cleaned/sPlot_plot_data.rds")
saveRDS(out, "data/cleaned/sPlot_trait_data.rds")
saveRDS(pres, "data/cleaned/sPlot_presence_data.rds")
saveRDS(trait2, "data/cleaned/trait_data_sPlot.rds")
