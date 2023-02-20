library(tidyverse)
library(Matrix)
library(doParallel)


d <- readRDS("data/cleaned/sPlot_trait_clim_data.rds") %>%
  mutate(plot_id = as.character(plot_id))


# function to generate presence-absence data

get_pres_abs <- function(df) {
  site_f <- factor(df$plot_id)
  species_f <- factor(df$species)
  
  m <- sparseMatrix(
    as.numeric(site_f), 
    as.numeric(species_f),
    x = 1,
    dimnames = list(levels(site_f), levels(species_f)))
  
  out <- tibble(plot_id = rownames(m)[row(m)],
                species = colnames(m)[col(m)],
                pres = as.numeric(m))
  
  return(out)
}



# presence-absence by ecoregion

d_eco <- split(d, d$ecoregion)

cl <- makeCluster(20)
registerDoParallel(cl)

pa_eco <- foreach(i = 1:length(d_eco), .combine = "rbind", .packages = c("dplyr", "Matrix")) %dopar% {
  get_pres_abs(d_eco[[i]])
}
stopCluster(cl)


# export data
saveRDS(pa_eco, "data/cleaned/pres_abs_eco.rds")
