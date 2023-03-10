library(tidyverse)
library(mgcv)

wtd_nt <- readRDS("results/models/arid_wtd_no_trait.rds")
wtd_wt <- readRDS("results/models/arid_wtd_traits.rds")

lrt <- anova(wtd_nt, wtd_wt, test = "Chisq")

aic_tab <- AIC(wtd_nt, wtd_wt) %>% 
  rownames_to_column() %>% 
  rename(model = 1) %>%
  arrange(AIC) %>%
  mutate(AICdif = AIC - AIC[1])

lrt <- lst(lrt, aic_tab)
saveRDS(lrt, "results/model_comparison.rds")
