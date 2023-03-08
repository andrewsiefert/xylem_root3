library(tidyverse)


d <- readRDS("data/cleaned/trait_data.rds") %>%
  as_tibble() %>%
  rename(species = species_wfo)


## TRY ----
try <- read_csv("~/projects/pheno/data/growthform/TRY/TRY_Categorical_Traits.csv") %>%
  janitor::clean_names() %>% 
  select(species = acc_species_name, pheno_try = leaf_phenology) %>%
  na.omit() %>%
  mutate(pheno_try = fct_recode(pheno_try, d = "deciduous", e = "evergreen", v = "deciduous/evergreen")) 

## Zanne ----

zanne <- read_csv("~/projects/pheno/data/growthform/Zanne/GlobalLeafPhenologyDatabase.csv") %>%
  rename(species = 1, pheno_z = 2) %>%
  mutate(pheno_z = fct_recode(pheno_z, d = "D", e = "EV", v = "D_EV")) 

## Combine ----
pheno <- left_join(d, try) %>%
  left_join(zanne) %>%
  mutate(leaf_phenology = case_when(
    pheno_try == "d" ~ "d",
    pheno_try == "e" ~ "e",
    pheno_z == "d" ~ "d",
    pheno_z == "e" ~ "e",
    pheno_try == "v" ~ "v", 
    pheno_z == "v" ~ "v", 
    TRUE ~ NA
  )) %>%
  select(species, leaf_phenology) %>%
  arrange(species)

write_csv(pheno, "data/traits/leaf_phenology.csv")
