library(tidyverse)
library(terra)


# read in aridity data
arid_r <- rast("data/clim/aridity/ai_et0.tif")

# read in sPlot plot data
plots <- readRDS("data/cleaned/sPlot_plot_data.rds")
coords <- distinct(plots, longitude, latitude)
points <- vect(as.matrix(coords), crs = crs(arid_r))

arid <- extract(arid_r, points)

out <- coords %>% 
  mutate(arid = arid$ai_et0/10000) %>%
  distinct() %>%
  na.omit()

saveRDS(out, "data/clim/aridity/sPlot_arid.rds")
