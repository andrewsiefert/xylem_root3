library(tidyverse)
library(ggExtra)


traits <- readRDS("data/cleaned/trait_data_sPlot.rds") 

d <- traits %>% 
  select(species, p50, rd_max) %>% 
  distinct() %>% 
  na.omit() 


# label angiosperms and gymnosperms
fams <- read_csv("data/phylo/phylo_species_list.csv")

d <- d %>% 
  left_join(fams) %>%
  mutate(p50_t = sqrt(-p50),
         clade = ifelse(family %in% c('Araucariaceae', 
                                      'Ginkgoaceae',
                                      'Pinaceae',
                                      'Sciadopityaceae',
                                      'Taxaceae',
                                      'Cupressaceae'),
                        'Gymnosperm', 
                        'Angiosperm'))

p50_labs <- -seq(0.5, 4, by = 0.5)^2
p50_breaks <- sqrt(-p50_labs)

# sample labeled species using bins and thresholds
ss <- d %>%
  sample_n(n()) %>%
  mutate(p50_bin = ntile(p50_t, 7), 
         rd_bin = ntile(log(rd_max), 7)) %>%
  group_by(p50_bin, rd_bin) %>%
  mutate(bin_n = 1:n()) %>%
  ungroup() %>%
  filter(bin_n == 1 | p50 > -1.5 | p50 < -6.5 | rd_max < 0.9 | rd_max > 12)

# quantiles
q <- expand.grid(p50_t = sqrt(c(1.5, 6.9)), 
                 rd_max = c(0.6, 12.9)) %>%
  mutate(label = c("Vulnerable confronters",
                   "Resistant confronters",
                   "Vulnerable avoiders",
                   "Resistant avoiders"), 
         clade = "Angiosperm")


p <- d %>% 
  ggplot(aes(x = p50_t, y = rd_max, shape = clade)) + 
  geom_hline(yintercept = median(d$rd_max), color = "gray") +
  geom_vline(xintercept = median(d$p50_t), color = "gray") +
  geom_point(alpha = 0.25, aes(color = clade)) + 
  scale_x_reverse(breaks = p50_breaks, labels = p50_labs) +
  scale_y_log10(limits = range(d$rd_max)) + 
  #geom_text(aes(label = species), fontface = 3, size = 3) +
  labs(y = "Rooting depth (m)") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = bquote(P[50]~(MPa)), shape = NULL, color = NULL) +
  ggrepel::geom_text_repel(data = ss, aes(label = species, color = clade), size = 3, max.overlaps = 5, 
                           fontface = 3, min.segment.length = 0.2, box.padding = 0.1, show.legend = F) +
  geom_point(data = ss, aes(color = clade)) +
  geom_point(data = q, size = 4, color = "gray40", show.legend = F) +
  scale_color_hue(l = 30) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12)) 

p


svg("results/figures/p50_rd_space.svg", height = 9.5, width = 9)
ggMarginal(p, size = 12, groupFill = T, color = NA)
dev.off()





# Older versions

# sample labeled species using bins
ss_d <- d %>%
  mutate(p50_bin = ntile(p50, 6), 
         rd_bin = ntile(log(rd_max), 6)) %>%
  group_by(p50_bin, rd_bin) %>%
  sample_n(1) %>%
  ungroup()

p <- d %>% 
  ggplot(aes(x = p50, y = rd_max, shape = clade, color = clade)) + 
  geom_hline(yintercept = median(d$rd_max), color = "gray") +
  geom_vline(xintercept = median(d$p50), color = "gray") +
  geom_point(alpha = 0.3) + 
  scale_x_continuous(limits = range(d$p50)) +
  scale_y_log10(limits = range(d$rd_max)) + 
  #geom_text(aes(label = species), fontface = 3, size = 3) +
  labs(y = "Rooting depth (m)") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "p50 (MPa)", shape = NULL, color = NULL) +
  ggrepel::geom_text_repel(data = ss_d, aes(label = species), size = 3, max.overlaps = 10, fontface = 3,
                           min.segment.length = 0.2, box.padding = 0.1, show.legend = F) +
  geom_point(data = ss_d) +
  scale_color_hue(l = 30) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12)) + 
  guides(shape = guide_legend(override.aes = list(size = 3))) 

#pdf("p50_rd_space_labels4.pdf", height = 9.5, width = 9)
ggMarginal(p, size = 10, fill = "gray", color = NA)
#dev.off()





dists <- d %>%
  mutate(rd = log(rd_max), 
         p50 = sqrt(abs(p50))) %>%
  select(p50, rd) %>%
  dist(diag = F, upper = T) %>% 
  as.matrix()
diag(dists) <- NA
dists <- apply(dists, 1, min, na.rm = T)

ss_d <- d %>%
  filter(dists > 0.24)

ss_b <- d %>%
  anti_join(ss_d)

bins <- anticlust::balanced_clustering(ss_b %>% select(p50, rd_max), 35)

ss_b <- ss_b %>%
  mutate(bin = bins) %>%
  filter(!is.na(bin)) %>%
  group_by(bin) %>%
  sample_n(1) %>%
  ungroup()

ss <- bind_rows(ss_b, ss_d) %>% distinct(species, p50, rd_max, clade)

p <- d %>% 
  ggplot(aes(x = p50, y = rd_max, shape = clade, color = clade)) + 
  geom_hline(yintercept = median(d$rd_max), color = "gray") +
  geom_vline(xintercept = median(d$p50), color = "gray") +
  geom_point(alpha = 0.3) + 
  scale_x_continuous(limits = range(d$p50)) +
  scale_y_log10(limits = range(d$rd_max)) + 
  #geom_text(aes(label = species), fontface = 3, size = 3) +
  labs(y = "Rooting depth (m)") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "p50 (MPa)", shape = NULL, color = NULL) +
  ggrepel::geom_text_repel(data = ss, aes(label = species), size = 3, max.overlaps = 5, fontface = 3,
                           min.segment.length = 0.2, box.padding = 0.1, show.legend = F) +
  geom_point(data = ss) +
  scale_color_hue(l = 30) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12)) + 
  guides(shape = guide_legend(override.aes = list(size = 3))) 

#pdf("p50_rd_space_labels5.pdf", height = 9.5, width = 9)
ggMarginal(p, size = 10, fill = "gray", color = NA)
#dev.off()


# older versions ----------------------------------------------------------

tr <- d %>% select(-species)


hihi <- hilo <- lohi <- lolo <- rep(NA, nrow(d))

for(i in 1:nrow(d)) {
  hihi[i] <- !any((d$p50[-i] >= d$p50[i]) & (d$rd_max[-i] >= d$rd_max[i])) 
  hilo[i] <- !any((d$p50[-i] >= d$p50[i]) & (d$rd_max[-i] <= d$rd_max[i])) 
  lohi[i] <- !any((d$p50[-i] <= d$p50[i]) & (d$rd_max[-i] >= d$rd_max[i])) 
  lolo[i] <- !any((d$p50[-i] <= d$p50[i]) & (d$rd_max[-i] <= d$rd_max[i])) 
}


tr_s <- tr %>% mutate(rd_max = log(rd_max)) %>% mutate_all(scale)

edge <- rowSums(cbind(hihi, hilo, lohi, lolo)) > 0
edge_spp <- d %>% filter(edge)

label <- rep(F, nrow(d))
label[order(as.matrix(dist(rbind(colMeans(tr_s), tr_s)))[-1,1])[1:3]] <- T
label[order(as.matrix(dist(rbind(colMeans(tr_s), tr_s)))[-1,1], decreasing = T)[1:6]] <- T
label[which.min(d$p50)] <- T
label[which.max(d$p50)] <- T

label_spp <- d %>% filter(label)


ggplot(d, aes(x = p50, y = rd_max)) + 
  geom_point(aes(color = label)) + 
  scale_y_log10() + 
  ggrepel::geom_label_repel(data = label_spp, aes(label = species)) +
  labs(y = "Rooting depth (m)") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "p50 (MPa)")

ggsave("p50_rd_space.pdf")


d %>% 
  ggplot(aes(x = p50, y = rd_max)) + 
  geom_point() + 
  scale_y_log10() + 
  ggrepel::geom_text_repel(aes(label = species), fontface = 3, size = 3, box.padding = 0.1, max.overlaps = 10,
                           min.segment.length = 0.1) +
  labs(y = "Rooting depth (m)") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "p50 (MPa)")
ggsave("p50_rd_space_labels2.pdf")    

d %>% 
  ggplot(aes(x = p50, y = rd_max)) + 
  geom_point(alpha = 0.2) + 
  scale_x_continuous(limits = range(d$p50)) +
  scale_y_log10(limits = range(d$rd_max)) + 
  #geom_text(aes(label = species), fontface = 3, size = 3) +
  labs(y = "Rooting depth (m)") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "p50 (MPa)") +
  ggrepel::geom_text_repel(data = ss, aes(label = species), size = 3, max.overlaps = 5, fontface = 3) +
  geom_point(data = ss)

ggsave("p50_rd_space_labels2.pdf")  

