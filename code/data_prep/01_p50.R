library(tidyverse)
library(WFO)
library(janitor)


# load Martin St.Paul et al. Ecol Lett 2017
martin <- read_csv("data/p50/martin.csv") %>%
  clean_names() %>% 
  select(species_binomial, group, p50, references) %>%
  rename(species = species_binomial, reference = references) %>%
  mutate(dataset = "martin") %>%
  filter(is.finite(p50)) %>%
  distinct()


# load new dataset compiled by Laughlin
new <- read_csv("data/p50/new.p50.data.filtered.csv") %>%
  clean_names() %>%
  rename(reference = paper) %>%
  select(species, group, p50, reference) %>%
  mutate(dataset = "laughlin") %>%
  filter(is.finite(p50)) %>%
  distinct()


# add XFT data
xft <- read_csv("data/p50/30MAR2021_cleaned_xft.csv") %>%
  as_tibble() %>%
  clean_names() %>%
  mutate(group = str_replace_all(group, c("gymnosperm"="Gymnosperm", "angiosperm"="Angiosperm"))) %>% 
  filter(plant_organ == "S", curve == "S", p50_method != "AD", p50_method != "AS") %>%
  select(species = binomial_gnr, group, p50, reference = short_reference_authors_year) %>%
  mutate(dataset = "xft") %>%
  filter(is.finite(p50)) %>%
  distinct()

xylem <- bind_rows(martin, new, xft) %>%
  filter(group %in% c("Angiosperm", "Gymnosperm"), 
         !str_detect(species, "sp\\.")) %>%
  arrange(species)


# standardize species names according to WFO taxonomic backbone

spp <- unique(xylem$species)

spp_list <- split(spp, ceiling(seq_along(spp)/50))
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
  select(species = spec.name.ORIG, species_wfo = scientificName) %>%
  as_tibble() %>% 
  full_join(distinct(xylem, species)) %>%
  mutate(flag = as.numeric(species != species_wfo))

write_csv(spp_lookup, "data/traits/p50_species_lookup.csv")
# make manual changes
spp_lookup <- read_csv("data/traits/p50_species_lookup.csv") %>% select(-flag)

xylem2 <- full_join(xylem, spp_lookup) %>%
  select(species, species_wfo, everything())

# calculate species average P50
p50 <- xylem2 %>% 
  group_by(species_wfo) %>%
  summarize(p50 = mean(p50)) %>%
  ungroup()

# export data
saveRDS(xylem2, "data/traits/p50_database.rds")
write_csv(p50, "data/cleaned/p50_data.csv")
