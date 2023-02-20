library(tidyverse)
library(V.PhyloMaker)
library(WorldFlora)
library(doParallel)


# Read in trait data ------------------------------------------------------

traits <- readRDS("data/cleaned/trait_data.rds") 

d <- traits %>% 
  select(species_wfo, p50, rd_max) %>% 
  distinct() %>% 
  na.omit() 



# Get standardized species names and families -----------------------------

species_list <- d$species_wfo


# read in taxonomic backbone
wfo <- read.delim("data/wfo_list.txt") %>% 
  filter(infraspecificEpithet == "")

# split up species list to allow parallelization
spp_list <- split(species_list, rep(1:11, each = 20))

# match species list to WFO taxonomic backbone
n <- length(spp_list)
cl <- makeCluster(n)
registerDoParallel(cl)

spp_match <- foreach(i = 1:n, .errorhandling = "pass", .packages = "WorldFlora") %dopar% {
  WFO.match(spp_list[[i]], WFO.data = wfo, Genus = "new_genus", Species = "new_species",
            exclude.infraspecific = T)
}

stopCluster(cl)

# create standardized species list
phylo_list <- spp_match %>% 
  map(WFO.one) %>%
  bind_rows() %>%
  select(species = scientificName, genus, family) %>%
  as_tibble() 

write_csv(phylo_list, "data/phylo/phylo_species_list.csv")

# Make the tree -----------------------------------------------------------

tree <- phylo.maker(sp.list = phylo_list, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S3")

write.tree(tree$scenario.3, "data/phylo/p50_rd_phylo.tre")
