library(tidyverse)


# read in sPlot plot data
plots <- readRDS("data/cleaned/sPlot_plot_data.rds")
d <- readRDS("data/cleaned/sPlot_trait_data.rds")

# read in climate data
arid <- read_rds("data/clim/aridity/sPlot_arid.rds")
wtd <- read_rds("data/clim/wtd/sPlot_wtd.rds")

clim <- inner_join(arid, wtd)

# merge in climate data
d2 <- d %>% 
  inner_join(plots) %>% 
  inner_join(clim) %>% 
  add_count(ecoregion) %>% 
  filter(n>=100) %>%
  select(-n)
plots2 <- inner_join(plots, clim)

# save output
saveRDS(d2, "data/cleaned/sPlot_trait_clim_data.rds")
saveRDS(plots2, "data/cleaned/plot_clim_data.rds")
