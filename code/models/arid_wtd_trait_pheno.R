library(tidyverse)
library(mgcv)

source("code/transformers.R")


d <- readRDS("data/cleaned/p50_rd_data_ecoregion.rds")


# add leaf phenology
pheno <- read_csv("data/traits/leaf_phenology.csv")
traits <- readRDS("data/cleaned/trait_data.rds") %>%
  select(species = species_wfo, p50, rd_max) 

pheno2 <- inner_join(pheno, traits) %>%
  mutate(decid = as.numeric(leaf_phenology == "d")) %>%
  select(p50:decid)

ecoregion <- distinct(d, ecoregion) %>%
  mutate(eco = factor(ecoregion))

species <- distinct(d, species) %>%
  mutate(spp = factor(species))

d <- d %>% 
  left_join(ecoregion) %>% 
  left_join(species) %>% 
  left_join(pheno2) %>%
  select(pres, arid, wtd, p50, rd_max, spp, eco, decid) %>%
  na.omit() %>%
  mutate(arid = transform(arid, "arid_sqrt"),
         arid2 = arid^2,
         wtd = transform(wtd, "wtd_log"),
         wtd2 = wtd^2,
         p50 = transform(p50, "P50_sqrt"),
         p502 = p50^2,
         rd = transform(rd_max, "rd_log"),
         rd2 = rd^2) %>%
  select(-rd_max)

start <- Sys.time()

m <- bam(pres ~ decid*(arid*wtd*p50*rd + arid2 + wtd2 + p502 + rd2) + s(spp, bs = "re") + s(eco, bs = "re"), 
         family = "binomial", data = d, discrete = T)

end <- Sys.time()
print(end-start)

saveRDS(m, "results/models/arid_wtd_traits_pheno.rds")
