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


df.full <- read.csv('data/df.csv') %>%
           select(-X) %>%
           mutate_if(is.character, factor) %>%
           mutate(roi = as.factor(roi),
                  rel_time = index - 30,
                  app_factor = as.factor(app),
                  dist_um = dist * 0.16,
                  dist_group = as.factor(case_when(dist_um <= 15 ~ 'I',
                                                   (dist_um > 30) & (dist_um <= 50)  ~ 'II',
                                                   (dist_um > 65) & (dist_um <= 200) ~ 'III',
                                                   .default = '0')))

##### EXPLORATORY ANALYSIS #####
df.summary <- df.full %>%
  filter(dist_group != '0') %>%
  group_by(lab_id, app_factor, dist_group) %>%
  summarise(n_cell = n_distinct(id), n_roi = n_distinct(roi),
            max = max(df))
  

ggplot(data = df.full %>% filter(dist_group != '0'), aes(x = df, fill = dist_group)) +
  geom_density(alpha = .5) +
  facet_wrap(~lab_id)


##### AMPLITUDE FILTERING #####
df.max <- df.full %>%
  filter(dist_group != '0', app_factor != '30') %>%
  select(df, lab_id, roi, app_factor, dist_group) %>%
  group_by(lab_id, roi, app_factor, dist_group) %>%
  summarise(max_df = max(df)) %>%
  distinct() %>%
  ungroup() %>%
  mutate(app = as.numeric(as.character(app_factor))) %>%
  filter(max_df > 0.001)

ggplot(data = df.max, aes(x = max_df, fill = lab_id)) +
  geom_density(alpha = .5) +
  facet_wrap(~app, nrow = 3, ncol = 3)

ggplot(data = df.max, aes(x = app, y = max_df)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(aes(color = dist_group, group = dist_group),
               fun = median,
               geom = 'line', linewidth = 0.75) +
  stat_summary(aes(color = dist_group, group = dist_group),
               fun = median,
               geom = 'point', size = 2) +
  stat_summary(aes(color = dist_group, group = dist_group),
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'errorbar', linewidth = 0.5, width=0.05) +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 2.5, 5, 10, 20, 60, 120)) +
  scale_y_continuous(breaks = seq(0,10,0.1)) +
  facet_wrap(~lab_id, nrow = 3) +
  theme_minimal_hgrid()



##### STAT 20 #####
df.20 <- df.full %>%
  filter(app_factor == '20',
         lab_id != 'psd')


ggplot(data = df.20, aes(x = rel_time, y = df)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  stat_summary(aes(color = dist_group, group = dist_group),
               fun = median,
               geom = 'line', linewidth = 0.75) +
  stat_summary(aes(color = dist_group, group = dist_group),
               fun = median,
               geom = 'point', size = 2) +
  stat_summary(aes(fill = dist_group, group = dist_group),
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .25) +
  facet_wrap(~lab_id)



ggplot(data = df.20 %>% filter(rel_time %in% c(0, 2, 5, 20, 21, 23)),
       aes(x = dist_um, y = df,
           color = lab_id, shape = lab_id, fill = lab_id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = .5) +
  geom_smooth(method = 'lm') +
  scale_y_continuous(limits = c(-0.25, 0.5)) +
  facet_wrap(~as.factor(rel_time)) +
  theme_classic()




##### DECAY 20 #####
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
