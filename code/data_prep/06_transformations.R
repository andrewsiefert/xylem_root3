library(tidyverse)


traits <- readRDS("data/cleaned/trait_data_sPlot.rds") 
plots <- readRDS("data/cleaned/plot_clim_data.rds")


# traits
traits_s <- traits %>%
  mutate(species = factor(species),
         p50_abs_log_scaled = scale(log(abs(p50))),
         p50_abs_sqrt_scaled = scale(sqrt(abs(p50))),
         rd_log_scaled = scale(log(rd_max))) %>%
  select(species, p50_abs_log_scaled, p50_abs_sqrt_scaled, rd_log_scaled)

square <- function(x) x^2

trait_trans <- sapply(traits_s[,-1], function(i) c(attr(i, "scaled:center"), attr(i, "scaled:scale"))) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  rename(trait = 1, center = 2, scale = 3)

write_csv(trait_trans, "data/cleaned/trait_transformations.csv")


# climate
plots_s <- plots %>% 
  mutate(plot_id = factor(plot_id),
         ecoregion = factor(ecoregion) %>% as.numeric() %>% as.factor(),
         arid_log_scaled = scale(log(arid)),
         arid_sqrt_scaled = scale(sqrt(arid)),
         wtd_log_scaled = scale(log(wtd+0.01))) %>%
  select(plot_id, ecoregion, arid_log_scaled, arid_sqrt_scaled, wtd_log_scaled)

plot_trans <- sapply(plots_s[,3:5], function(i) c(attr(i, "scaled:center"), attr(i, "scaled:scale"))) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  rename(variable = 1, center = 2, scale = 3)

write_csv(plot_trans, "data/cleaned/plot_transformations.csv")
