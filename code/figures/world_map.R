library(tidyverse)
library(sf)

sf::sf_use_s2(FALSE)


dunia <- spData::world %>% 
  st_as_sf()

d <- readRDS("data/cleaned/sPlot_trait_clim_data.rds") %>%
  distinct(plot_id, longitude, latitude)

plots <- st_as_sf(d, coords = c("longitude","latitude"), crs = st_crs(dunia)) 


grid <- st_make_grid(plots, square = F, cellsize = 3)

tab <- st_intersects(grid, plots)

grid <- st_sf(n = lengths(tab), geometry = st_cast(grid, "MULTIPOLYGON")) 
grid <- grid %>% filter(n>0)

ggplot() + 
  geom_sf(data = dunia %>% st_transform("ESRI:54012"), size = 0.3) + 
  geom_sf(data = grid %>% st_transform("ESRI:54012"), aes(fill = n), alpha = 0.7, color = NA) +
  scale_fill_viridis_c(trans = "log10", breaks = c(10, 100, 1000, 10000)) + 
  cowplot::theme_map() +
  theme(panel.grid = element_line(color = "gray80")) +
  labs(fill = "# plots")

ggsave("results/figures/plot_map.jpg", height = 4, width = 7, units = "in", dpi = 300)
ggsave("results/figures/plot_map.svg", height = 4, width = 7, units = "in")
