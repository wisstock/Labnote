require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)
library(gridExtra)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_06_19_Sasa_HPCA-PSD')


df.full <- read.csv('data/all_cells_hpca_psd95.csv') %>%
           select(id, lab_id, roi, dist, index, time, abs_int, dF_int, df, app) %>%
           mutate_if(is.character, factor) %>%
           mutate(roi = as.factor(roi))

write.csv(df.full, "df.csv")


##### DECAY #####
tail.20 <- seq(50,65)
df.lm.20 <- df.full %>%
  filter(app == 20, time %in% tail.20, lab_id != 'psd') %>%
  select(time, df, lab_id, roi, id) %>%
  mutate(idx = row_number(),
         df = if_else(df < 0, 0.000000000000001, df)) %>%
  distinct() %>%
  droplevels() %>%
  select(roi, lab_id, id, df, idx) %>%
  group_by(roi, lab_id, id) %>%
  nest() %>%
  mutate(model = map(data, ~ tidy(lm(log(df) ~ idx, data = .x)))) %>%
  unnest(model) %>%
  ungroup() %>%
  select(-std.error, -statistic, -p.value) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  mutate(decay_tau = -1/idx)

  
ggplot(data = df.lm.20 %>% filter(decay_tau < 100, decay_tau > 1),
       aes(fill = lab_id, x = decay_tau)) +
  geom_density(alpha = .5)
