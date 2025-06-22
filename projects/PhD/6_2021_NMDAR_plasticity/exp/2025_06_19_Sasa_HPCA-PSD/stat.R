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
                  app_factor = as.factor(app),
                  dist_um = dist * 0.16,
                  dist_group = as.factor(case_when(dist_um <= 25 ~ 'I',
                                                   (dist_um > 25) & (dist_um <= 50)  ~ 'II',
                                                   dist_um > 50 ~ 'III',
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



##### AMPLITUDE 20 #####
df.20 <- df.full %>%
  filter(app_factor == '20') %>%
  mutate(rel_time = index - 30)


ggplot(data = df.20,
       aes(x = rel_time, y = df)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  stat_summary(aes(color = dist_group, group = dist_group),
               fun = median,
               geom = 'line', linewidth = 0.75) +
  stat_summary(aes(color = dist_group, group = dist_group),
               fun = median,
               geom = 'point', size = 1) +
  stat_summary(aes(fill = dist_group, group = dist_group),
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .25) +
  facet_wrap(~lab_id)

# cor test
df.20.cor <- df.20 %>%
  filter(rel_time %in% seq(-5,30)) %>%
  group_by(lab_id, rel_time, dist_group) %>%
  cor_test(dist_um, df, method = 'pearson')
  
ggplot(data = df.20.cor, aes(x = rel_time, y = cor, color = lab_id, fill = lab_id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = c(0, 20), linetype = 3) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              size = 0, alpha = .3) +
  labs(title = 'Pearson correlation (distance vs. amplitude) with 95% CI',
       x = 'Time (s)',
       y = 'Cor. coef.',
       color = 'Region',
       fill = 'Region') +
  facet_wrap(~dist_group) +
  theme_minimal() +
  theme(legend.position = c(0.1, 0.8))

# slope test
df.20.lm <- df.20 %>%
  filter(rel_time %in% seq(-5,30)) %>%
  select(lab_id, rel_time, df, dist_um, dist_group) %>%
  group_by(lab_id, rel_time, dist_group) %>%
  nest() %>%
  mutate(model = map(data, ~ tidy(lm(df ~ dist_um, data = .x)))) %>%
  unnest(model) %>%
  ungroup() %>%
  select(-data) %>%
  filter(term == 'dist_um') %>%
  mutate(signig = as.factor(case_when(p.value < 0.05 ~ '*',
                                      p.value > 0.05 ~ 'ns')))

ggplot(data = df.20.lm, aes(x = rel_time, y = estimate, color = lab_id, fill = lab_id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = c(0, 20), linetype = 3) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = estimate-(std.error*1.96), ymax = estimate+(std.error*1.96)),
              size = 0, alpha = .3) +
  facet_wrap(~dist_group) +
  labs(title = 'Linear fit slope with 95% CI',
       x = 'Time (s)',
       y = 'Slope',
       color = 'Region',
       fill = 'Region') +
  theme_minimal() +
  theme(legend.position = c(0.1, 0.8))
  

# selected frames lm
time.points.20 <- c(0, 1, 2, 3, 4, 5)

ggplot(data = df.20 %>%
              filter(rel_time %in% time.points.20, lab_id != 'shaft') %>%
              mutate(lab_dist = interaction(lab_id, dist_group, sep = '_')),
       aes(x = dist_um, y = df,
           color = lab_id, shape = lab_id, fill = lab_id, group = lab_dist)) +
  geom_vline(xintercept = c(25, 50), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = .1) +
  geom_smooth(method = 'lm') +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  facet_wrap(~as.factor(rel_time)) +
  labs(title = 'Linear fit for selected frames (app. 20 s)',
       x = 'Distance (μm)',
       y = 'ΔF/F0',
       color = 'Region',
       fill = 'Region',
       shape = 'Region') +
  theme_minimal()


##### 20 BOX STAT #####


##### 20 LINE STAT #####

df.20.for.stat <- df.20 %>%
  mutate(rel_time = index - 30) %>%
  filter(rel_time %in% seq(-2,25)) %>%
  select(df, lab_id, dist_group, id, rel_time)

df.20.stat <- df.20.for.stat %>%
  filter(lab_id != 'psd') %>%
  group_by(dist_group, rel_time) %>%
  wilcox_test(df ~ lab_id) %>%
  add_significance() %>%
  add_y_position(fun = 'median')

ggplot(data = df.20.for.stat,
       aes(x = rel_time, y = df)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(aes(color = lab_id, group = lab_id),
               fun = median,
               geom = 'line', linewidth = 0.75) +
  stat_summary(aes(color = lab_id, group = lab_id),
               fun = median,
               geom = 'point', size = 1) +
  stat_summary(aes(fill = lab_id, group = lab_id),
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .25) +
  geom_text(data = df.20.stat %>% filter(p.signif != 'ns'),
            aes(label = p.signif, x = rel_time, y = y.position-0.02)) +
  facet_wrap(~dist_group)
  facet_wrap(~interaction(dist_group, rel_time, sep = '_'), ncol = 3) +
  theme_minimal() +
  theme(legend.position = 'none')
  
  

##### DECAY 20 #####
tail.20 <- seq(20,30)
df.lm.20 <- df.20 %>%
  filter(rel_time %in% tail.20, lab_id != 'psd') %>%
  select(time, df, lab_id, roi, id, dist_group) %>%
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
