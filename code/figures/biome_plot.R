library(tidyverse)
library(plotbiomes)
library(terra)
library(prioritizr)
library(doParallel)


# read in data
d <- readRDS("data/cleaned/sPlot_presence_data.rds") %>% as_tibble()
traits <- readRDS("data/cleaned/trait_data_sPlot.rds") 
plots <- readRDS("data/cleaned/plot_clim_data.rds") %>% dplyr::select(plot_id:ecoregion)


keep_plots <- d %>% 
  inner_join(traits) %>% 
  filter(is.finite(p50), is.finite(rd_max)) %>%
  distinct(plot_id)

# merge data
d2 <- plots %>% semi_join(keep_plots)

rm(d)
rm(plots)
rm(traits)


# extract MAT and annual precip
mat_r <- rast("data/clim/worldclim/wc2.1_2.5m_bio_1.tif")
ppt_r <- rast("data/clim/worldclim/wc2.1_2.5m_bio_12.tif")

coords <- distinct(d2, longitude, latitude)

mat <- extract(mat_r, coords)
ppt <- extract(ppt_r, coords)

coords$mat <- mat[,2]
coords$ppt <- ppt[,2]

d3 <- left_join(d2, coords) %>% 
  mutate(ppt = round(ppt/10)) %>% 
  filter(mat > -20, ppt < 650) %>% 
  mutate(mat_bin = cut_interval(mat, n = 40),
         ppt_bin = cut_interval(ppt, n = 40)) %>%
  group_by(mat_bin) %>%
  mutate(mat = (max(mat)+min(mat))/2) %>%
  ungroup() %>%
  group_by(ppt_bin) %>%
  mutate(ppt = (max(ppt)+min(ppt))/2) %>%
  ungroup() %>%
  group_by(mat_bin, ppt_bin) %>%
  summarize(mat = mean(mat), 
            ppt = mean(ppt), 
            n = n()) %>%
  ungroup()


whittaker <- Whittaker_biomes %>% rename(mat = 1, ppt = 2)


ggplot(whittaker, aes(x = mat, y = ppt)) + 
  geom_polygon(aes(fill = biome)) +
  geom_point(data = d3, aes(size = n), alpha = 0.5) +
  scale_fill_manual(name   = "Whittaker biomes",
                    breaks = names(Ricklefs_colors),
                    labels = names(Ricklefs_colors),
                    values = Ricklefs_colors) +
  scale_size(trans = "sqrt", breaks = c(1, 100, 10000), range = c(0.5, 3)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = bquote(Temperature~(degree*C)), y = "Precipitation (cm)", size = "# plots")

ggsave("results/figures/biome_plot.svg", height = 5, width = 7.5, units = "in")
