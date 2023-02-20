library(tidyverse)
library(mgcv)
library(terra)
library(ks)
library(modelr)


source("code/transformers.R")



# Read in model and data --------------------------------------------------

m <- readRDS("results/models/arid_wtd_traits.rds")

d <- readRDS("data/cleaned/p50_rd_data_ecoregion.rds")

ecoregion <- distinct(d, ecoregion) %>%
  mutate(eco = factor(ecoregion))

species <- distinct(d, species) %>%
  mutate(spp = factor(species))

d <- d %>% 
  left_join(ecoregion) %>% 
  left_join(species) %>% 
  select(pres, arid, wtd, p50, rd_max, spp, eco) %>%
  na.omit() %>%
  mutate(arid = transform(arid, "arid_sqrt"),
         arid2 = arid^2,
         wtd = transform(wtd, "wtd_log"),
         wtd2 = wtd^2,
         p50 = transform(p50, "P50_sqrt"),
         p502 = p50^2,
         rd = transform(rd_max, "rd_log"),
         rd2 = rd^2) %>%
  dplyr::select(-rd_max)

trait <- readRDS("data/cleaned/trait_data_sPlot.rds") %>%
  dplyr::select(species, p50, rd_max) %>%
  na.omit()

plot <- readRDS("data/cleaned/plot_clim_data.rds")

# Prepare stuff for plotting ----------------------------------------------

# get posterior samples
post <- rmvn(2000, coef(m)[1:20], vcov(m)[1:20,1:20]) %>%
  as_tibble() %>%
  setNames(names(coef(m))[1:20]) %>%
  janitor::clean_names()

# calculate intercepts accounting for weighted mean random effects
coef <- coef(m)
eco <- coef[str_detect(names(coef), "eco")]
spp <- coef[str_detect(names(coef), "spp")]

eco_int <- weighted.mean(eco, table(d$eco))
spp_int <- weighted.mean(spp, table(d$spp))
re_int <- eco_int + spp_int

# calculate trait & environment sequences and quantiles
arid_r <- quantile(d$arid, c(0.01, 0.99))
arid_seq <- seq(arid_r[1], arid_r[2], length.out = 100)
arid_levs <- quantile(plot$arid, c(0.05, 0.95)) %>% transform("arid_sqrt")

wtd_r <- quantile(d$wtd, c(0.01, 0.99))
wtd_seq <- seq(wtd_r[1], wtd_r[2], length.out = 100)
wtd_levs <- quantile(plot$wtd, c(0.05, 0.95)) %>% transform("wtd_log")

p50_seq <- quantile(d$p50, seq(0.01, 0.99, length.out = 100))
p50_levs <- quantile(trait$p50, c(0.05, 0.95)) %>% transform("P50_sqrt")

p50_labs <- c(-2, -4, -6, -8)
p50_breaks <- transform(p50_labs, "P50_sqrt")

rd_seq <- quantile(d$rd, seq(0.01, 0.99, length.out = 100))
rd_levs <- quantile(trait$rd_max, c(0.05, 0.95)) %>% transform("rd_log")


# axis labels
p50_lab <- bquote(P[50]~(MPa))
rd_lab <- "Rooting depth (m)"
wtd_lab <- "Water table depth (m)"


# labeling functions
label_p50 <- function(x, bt = TRUE) {
  if(bt) x <- backtransform(x, "P50_sqrt")
  lab <- x %>% round(1) %>% paste("P50 =", ., "MPa") %>% fct_reorder(x)
  return(lab)
}

label_rd <- function(x, bt = TRUE) {
  if(bt) x <- backtransform(x, "rd_log")
  lab <- x %>% round(1) %>% paste("Rooting depth =", ., "m") %>% fct_reorder(-x)
  return(lab)
}

label_arid <- function(x, bt = TRUE) {
  if(bt) x <- backtransform(x, "arid_sqrt")
  lab <- x %>% round(1) %>% paste("Aridity Index =", .) %>% fct_reorder(x)
  return(lab)
}

label_wtd <- function(x, bt = TRUE) {
  if(bt) x <- backtransform(x, "wtd_log")
  lab <- x %>% round(1) %>% paste("WTD =", ., "m") %>% fct_reorder(-x)
  return(lab)
}


# colors for masks in 2D plots

t_col <- function(color, percent = 50) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, 
               alpha = (100 - percent) * 255 / 100)
  invisible(t.col)
}

clear <- t_col("white", 100)
mid <- t_col("white", 25)


# Single trait effects -----------------------------------------------------------

## P50 ----

grid <- expand.grid(p50 = p50_seq,
                    arid = arid_levs,
                    wtd = wtd_levs) %>%
  as_tibble() %>%
  mutate(p502 = p50^2, 
         y = 1)

X <- model.matrix(y ~ arid*wtd*p50 + p502, data = grid) %>% 
  as_tibble() %>% 
  janitor::clean_names() %>%
  as.matrix()

beta <- t(post[colnames(X)])

pred <- plogis(X %*% beta + re_int)

df <- grid %>%
  mutate(arid = label_arid(arid),
         wtd = label_wtd(wtd),
         est = apply(pred, 1, mean),
         lo = apply(pred, 1, quantile, 0.05),
         hi = apply(pred, 1, quantile, 0.95))

labs <- df %>%
  mutate(p50 = max(p50), est = max(hi)*1.05) %>%
  distinct(arid, wtd, p50, est) %>%
  mutate(panel = c("C", "D", "A", "B"))

ggplot(df, aes(x = p50, y = est)) + 
  geom_path() + 
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_text(data = labs, aes(label = panel), fontface = "bold", size = 4) +
  scale_x_reverse(breaks = p50_breaks, labels = p50_labs) +
  scale_y_continuous(expand = expansion(mult = 0.07)) +
  facet_grid(wtd~arid) +
  labs(x = p50_lab, 
       y = "Occurrence probability") +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave("results/figures/p50_effect_arid_wtd.jpg", height = 4.5, width = 6, units = "in")


## Rooting depth ----

grid <- expand.grid(rd = rd_seq,
                    arid = arid_levs,
                    wtd = wtd_levs) %>%
  as_tibble() %>%
  mutate(rd2 = rd^2, 
         y = 1)

X <- model.matrix(y ~ arid*wtd*rd + rd2, data = grid) %>% 
  as_tibble() %>% 
  janitor::clean_names() %>%
  as.matrix()

beta <- t(post[colnames(X)])

pred <- plogis(X %*% beta + re_int)

df <- grid %>%
  mutate(rd = backtransform(rd, "rd_log"),
         arid = label_arid(arid),
         wtd = label_wtd(wtd),
         est = apply(pred, 1, mean),
         lo = apply(pred, 1, quantile, 0.05),
         hi = apply(pred, 1, quantile, 0.95))

labs <- df %>%
  mutate(rd = min(rd), est = max(hi)*1.05) %>%
  distinct(arid, wtd, rd, est) %>%
  mutate(panel = c("C", "D", "A", "B"))

ggplot(df, aes(x = rd, y = est)) + 
  geom_path() + 
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_text(data = labs, aes(label = panel), fontface = "bold", size = 4.5) +
  facet_grid(wtd~arid) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(mult = 0.07)) +
  labs(x = rd_lab, 
       y = "Occurrence probability") +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave("results/figures/rd_effect_arid_wtd.jpg", height = 4.5, width = 6, units = "in")


# Single environment responses ----------------------------------------------

## Aridity response ----

grid <- expand.grid(p50 = p50_levs,
                    rd = rd_levs, 
                    arid = arid_seq) %>%
  as_tibble() %>%
  mutate(arid2 = arid^2, 
         y = 1)

X <- model.matrix(y ~ arid*p50*rd + arid2, data = grid) %>% 
  as_tibble() %>% 
  janitor::clean_names() %>%
  as.matrix()

beta <- t(post[colnames(X)])

pred <- plogis(X %*% beta + re_int)

df <- grid %>%
  mutate(p50_b = label_p50(p50),
         rd_b = label_rd(rd),
         arid = backtransform(arid, "arid_sqrt"),
         est = apply(pred, 1, mean),
         lo = apply(pred, 1, quantile, 0.05),
         hi = apply(pred, 1, quantile, 0.95))

labs <- df %>%
  mutate(arid = min(arid), est = max(hi)) %>%
  distinct(arid, p50_b, rd_b, est) %>%
  mutate(panel = c("C", "D", "A", "B"))

ggplot(df, aes(x = arid, y = est)) + 
  geom_path() + 
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_text(data = labs, aes(label = panel), fontface = "bold", size = 4) +
  scale_y_continuous(expand = expansion(mult = 0.07)) +
  facet_grid(rd_b~p50_b) +
  labs(x = "Aridity Index", 
       y = "Occurrence probability") +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave("results/figures/aridity_response_arid_wtd.jpg", height = 4.5, width = 6, units = "in")


## WTD response ---- 

grid <- expand.grid(p50 = p50_levs,
                    rd = rd_levs, 
                    wtd = wtd_seq) %>%
  as_tibble() %>%
  mutate(wtd2 = wtd^2, 
         y = 1)

X <- model.matrix(y ~ wtd*p50*rd + wtd2, data = grid) %>% 
  as_tibble() %>% 
  janitor::clean_names() %>%
  as.matrix()

beta <- t(post[colnames(X)])

pred <- plogis(X %*% beta + re_int)

df <- grid %>%
  mutate(p50_b = label_p50(p50),
         rd_b = label_rd(rd),
         wtd = backtransform(wtd, "wtd_log")+0.01,
         est = apply(pred, 1, mean),
         lo = apply(pred, 1, quantile, 0.05),
         hi = apply(pred, 1, quantile, 0.95))

labs <- df %>%
  mutate(wtd = min(wtd), est = max(hi)) %>%
  distinct(wtd, p50_b, rd_b, est) %>%
  mutate(panel = c("C", "D", "A", "B"))

ggplot(df, aes(x = wtd, y = est)) + 
  geom_path() + 
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_text(data = labs, aes(label = panel), fontface = "bold", size = 4) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(mult = 0.07)) +
  facet_grid(rd_b~p50_b) +
  labs(x = wtd_lab, 
       y = "Occurrence probability") +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave("results/figures/wtd_response_arid_wtd.jpg", height = 4.5, width = 6, units = "in")


# 2D environmental response ---------------------------------------------------

env <- distinct(d, arid, wtd)

n <- 500
bw <- Hpi(cbind(env$arid, env$wtd))*4

de <- kde(cbind(env$arid, env$wtd), H = bw, gridsize = n)

r <- terra::rast(t(de$estimate)[n:1,])
ext(r) <- c(min(env$arid), max(env$arid), min(env$wtd), max(env$wtd))

cont <- as.contour(r, nlevels = 100)

poly <- as.polygons(cont[cont$level == cont$level[2]])
poly2 <- as.polygons(cont[cont$level == cont$level[1]])

r2 <- rast(ncols=200, nrows=200, extent = ext(poly2), vals = 1)
mask <- rasterize(poly, r2, fun = sum) %>%
  as.data.frame(xy = T, na.rm = F) %>%
  as_tibble() %>%
  rename(arid = 1, wtd = 2, mask = 3)

grid <- expand.grid(p50 = p50_levs,
                    rd = rd_levs,
                    arid = seq_range(mask$arid, 50),
                    wtd = seq_range(mask$wtd, 50)) %>%
  na.omit() %>%
  #filter(arid < 3.5) %>%
  as_tibble() %>%
  mutate(arid2 = arid^2,
         wtd2 = wtd^2, 
         p502 = p50^2,
         rd2 = rd^2,
         y = 1)

X <- model.matrix(y ~ arid*wtd*p50*rd + arid2 + wtd2 + p502 + rd2, data = grid) %>% 
  as_tibble() %>% 
  janitor::clean_names() %>%
  as.matrix()

beta <- t(post[colnames(X)])

pred <- plogis(X %*% beta + re_int) %>% rowMeans()

df <- grid %>%
  mutate(p50_b = label_p50(p50),
         rd_b = label_rd(rd),
         arid = backtransform(arid, "arid_sqrt"),
         wtd = backtransform(wtd, "wtd_log") + 0.01,
         pres = pred)

mask <- mask %>%
  mutate(arid = backtransform(arid, "arid_sqrt"),
         wtd = backtransform(wtd, "wtd_log") + 0.01)

labs <- df %>%
  mutate(arid = min(arid), wtd = max(wtd)) %>%
  distinct(arid, wtd, rd_b, p50_b) %>%
  mutate(panel = c("C", "D", "A", "B"))

ggplot(df, aes(x = arid, y = wtd)) + 
  geom_raster(aes(fill = pres)) +
  geom_contour(aes(z = pres), color = "white", breaks = seq(0, 1, 0.01)) +
  scale_fill_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1)) +  
  labs(x = "Aridity index", 
       y = "Water table depth (m)",
       fill = "Probability of\noccurrence") +
  ggnewscale::new_scale_fill() +
  geom_raster(data = mask, aes(x = arid, y = wtd, fill = mask), show.legend = F) +
  geom_text(data = labs, aes(label = panel), fontface = "bold", size = 4.5, 
            nudge_x = 0.07, nudge_y = -0.15) +
  scale_fill_gradient2(low = clear, mid = clear, high = clear, na.value = mid) + 
  facet_grid(rd_b~p50_b) +
  scale_x_sqrt(expand = c(0,0), limits = range(df$arid), oob = scales::squish_infinite) +
  scale_y_log10(expand = c(0,0), limits = range(df$wtd), oob = scales::squish_infinite,) +
  theme_bw()

ggsave("results/figures/env2d_arid_wtd.jpg", height = 6, width = 7.5, dpi = 600)
ggsave("results/figures/env2d_arid_wtd.svg", height = 6, width = 7.5, units = "in")

# 2D trait effect ---------------------------------------------------------

# create mask

tr <- d %>% 
  dplyr::select(p50, rd) %>% 
  distinct() %>%
  mutate(p50 = sqrt(abs(backtransform(p50, "P50_sqrt"))),
         rd = log(backtransform(rd, "rd_log")))

hull_pts <- chull(tr)
hull <- tr[c(hull_pts, hull_pts[1]),] %>% 
  as_tibble()

p <- sp::Polygon(hull)
ps <- sp::Polygons(list(p),1)
sps <- sp::SpatialPolygons(list(ps))


mask <- expand.grid(x = modelr::seq_range(tr$p50, 300), 
                    y = modelr::seq_range(tr$rd, 300)) %>%
  mutate(value = 1) %>%
  raster::rasterFromXYZ() %>%
  raster::crop(raster::extent(sps)) %>%
  raster::mask(sps, inverse = T) %>%
  as.data.frame(xy = T) %>%
  as_tibble() %>%
  mutate(x = x, y = exp(y))

# get predictions
grid <- expand.grid(arid = arid_levs, 
                    wtd = wtd_levs,
                    p50 = modelr::seq_range(d$p50, 50, expand = 0.01), 
                    rd = modelr::seq_range(transform(exp(tr$rd), "rd_log"), 50, expand = 0.01)) %>%
  as_tibble() %>%
  mutate(arid2 = arid^2,
         wtd2 = wtd^2,
         p502 = p50^2,
         rd2 = rd^2, 
         y = 1)

X <- model.matrix(y ~ arid*wtd*p50*rd + arid2 + wtd2 + p502 + rd2, data = grid) %>% 
  as_tibble() %>% 
  janitor::clean_names() %>%
  as.matrix()

beta <- t(post[colnames(X)])

pred <- plogis(X %*% beta + re_int) %>% rowMeans()

df <- grid %>%
  mutate(p50 = backtransform(p50, "P50_sqrt") %>% abs() %>% sqrt(),
         rd = backtransform(rd, "rd_log"),
         arid = label_arid(arid),
         wtd = label_wtd(wtd),
         pres = pred)

labs <- df %>%
  mutate(p50 = max(p50), rd = max(rd)) %>%
  distinct(arid, wtd, rd, p50) %>%
  mutate(panel = c("C", "D", "A", "B"))

p50_labs <- -seq(0.5, 4, by = 0.5)^2
p50_breaks <- sqrt(-p50_labs)

df %>% 
  ggplot(aes(x = p50, y = rd)) + 
  geom_raster(aes(fill = pres %>% pmax(0.0001) %>% pmin(0.3))) + 
  geom_contour(aes(z = pres), color = "white", breaks = seq(0, 1, 0.02)) +
  scale_fill_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1)) +   
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = p50_lab, y = rd_lab, fill = "Probability of\noccurrence") +
  ggnewscale::new_scale_fill() +
  geom_raster(data = mask, aes(x = x, y = y, fill = value), show.legend = F) +
  geom_text(data = labs, aes(label = panel), fontface = "bold", size = 4.5, 
            nudge_x = -0.15, nudge_y = -0.11) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(wtd~arid) +
  geom_point(data = tr %>% mutate(rd = exp(rd)), alpha = 0.2) +
  scale_x_continuous(expand=c(0,0), limits = c(max(tr$p50), min(tr$p50)),
                  oob = scales::squish_infinite, labels = p50_labs, breaks = p50_breaks) +
  scale_y_log10(expand=c(0,0), limits = range(df$rd), oob = scales::squish_infinite) 

ggsave("results/figures/trait2d_arid_wtd.jpg", height = 6, width = 7.5, units = "in", dpi = 600)
ggsave("results/figures/trait2d_arid_wtd.svg", height = 6, width = 7.5, units = "in")
