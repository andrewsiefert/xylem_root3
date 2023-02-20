library(tidyverse)


# read in data
d <- readRDS("data/cleaned/pres_abs_eco.rds") %>% as_tibble()
traits <- readRDS("data/cleaned/trait_data_sPlot.rds") 
plots <- readRDS("data/cleaned/plot_clim_data.rds") 


# merge data
d <- d %>%
  mutate(plot_id = as.numeric(plot_id)) %>%
  inner_join(traits) %>%
  inner_join(plots)

rm(list = setdiff(ls(), "d"))


# P50 + rooting depth data ------------------------------------------------

out <- d %>%
  select(species, ecoregion, p50, rd_max, arid, wtd, pres) %>%
  mutate(species = as.numeric(factor(species)),
         ecoregion = as.numeric(factor(ecoregion))) 

saveRDS(out, "data/cleaned/p50_rd_data_ecoregion.rds")
