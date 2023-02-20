library(tidyverse)
library(terra)


# load and merge WTD rasters
paths <- list.files("data/clim/wtd", full.names = T) %>% str_subset("annualmean")

r <- map(paths, rast) %>%
  map(~mask(.x$WTD, .x$mask, maskvalues = 0))

wtd_r <- merge(r[[1]], r[[2]], r[[3]], r[[4]], r[[5]])
rm(r)


# read in sPlot plot data
plots <- readRDS("data/cleaned/sPlot_plot_data.rds")
coords <- distinct(plots, longitude, latitude)

wtd <- extract(wtd_r, coords) %>% mutate(WTD = -WTD/10)

out <- coords %>% 
  mutate(wtd = wtd$WTD) %>% 
  distinct() %>%
  na.omit()

saveRDS(out, "data/clim/wtd/sPlot_wtd.rds")
