library(tidyverse)
library(WorldFlora)
library(doParallel)


# load data
rooting_depth <- read_csv("data/traits/rsip_tumbar.csv") %>% rename(species = Species)

rd <- rooting_depth %>% 
  rename_all(tolower) %>%
  mutate(binomial = str_extract(species, "[A-Za-z]+ [A-Za-z-]+")) %>%
  select(binomial, growth_form, dr) %>%
  arrange(binomial) %>%
  mutate(growth_form = tolower(growth_form)) %>%
  filter(dr > 0.1, 
         growth_form %in% c("tree", "shrub", "semi-shrub"),
         !str_detect(binomial, "spp|sp$|sp\\.|^[a-z]"))

# get accepted species names
spp <- unique(rd$binomial)

spp_list <- split(spp, ceiling(seq_along(spp)/62))
length(spp_list)

wfo <- read.delim("data/wfo_list.txt")

n <- length(spp_list)
cl <- makeCluster(n)
registerDoParallel(cl)

spp_std <- foreach(i = 1:n, .errorhandling = "pass", .packages = "WorldFlora") %dopar% {
  WFO.match(spp_list[[i]], WFO.data = wfo, exclude.infraspecific = T)
}
stopCluster(cl)

spp_lookup <- spp_std %>%
  map(WFO.one) %>%
  bind_rows() %>%
  select(binomial = spec.name.ORIG, species_wfo = scientificName) %>%
  as_tibble() %>% 
  full_join(distinct(rd, binomial)) %>%
  mutate(flag = as.numeric(binomial != species_wfo))

write_csv(spp_lookup, "data/traits/rd_species_lookup.csv")
# make manual changes
spp_lookup <- read_csv("data/traits/rd_species_lookup.csv") %>% select(-flag)

# Merge in accepted names and collapse duplicates ----

rd2 <- rd %>%
  inner_join(spp_lookup) %>%
  select(species_wfo, dr)

rd_means <- rd2 %>%
  group_by(species_wfo) %>%
  summarize(rd_mean = mean(dr), rd_max = max(dr)) %>%
  ungroup()

# Merge data for woody species----

p50 <- read_csv("data/cleaned/p50_data.csv")

traits <- rd_means %>% 
  inner_join(p50) %>%
  select(species_wfo, p50, rd_mean, rd_max)


# export data
saveRDS(traits, "data/cleaned/trait_data.rds")
saveRDS(rd2, "data/cleaned/rooting_data.rds")

