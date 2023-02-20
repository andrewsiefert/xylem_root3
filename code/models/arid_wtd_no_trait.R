library(tidyverse)
library(mgcv)

source("code/transformers.R")


d <- readRDS("data/cleaned/p50_rd_data_ecoregion.rds")

ecoregion <- distinct(d, ecoregion) %>%
  mutate(eco = factor(ecoregion))

species <- distinct(d, species) %>%
  mutate(spp = factor(species))

d <- d %>% 
  left_join(ecoregion) %>% 
  left_join(species) %>% 
  select(pres, arid, wtd, spp, eco) %>%
  na.omit() %>%
  mutate(arid = transform(arid, "arid_sqrt"),
         wtd = transform(wtd, "wtd_log"), 
         arid2 = arid^2,
         wtd2 = wtd^2) 

start <- Sys.time()

m <- bam(pres ~ arid*wtd + arid2 + wtd2 + s(spp, bs = "re") + s(eco, bs = "re"), 
         family = "binomial", data = d,
         discrete = T)

end <- Sys.time()
print(end-start)

saveRDS(m, "results/models/arid_wtd_no_trait.rds")
