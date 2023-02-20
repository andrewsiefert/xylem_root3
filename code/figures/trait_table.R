library(tidyverse)


# read in species mean traits
traits <- read_csv("data/cleaned/trait_data_sPlot.csv")


# add observation counts
p50_n <- readRDS("data/traits/p50_database.rds") %>%
  count(species_wfo) %>%
  rename(species = 1, n_p50 = n)

rd_n <- readRDS("data/cleaned/rooting_data.rds") %>%
  count(species_wfo) %>%
  rename(species = 1, n_rd = 2)

out <- left_join(traits, p50_n) %>% left_join(rd_n)

write_csv(out, "results/tables/trait_table.csv")
