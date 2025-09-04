# Sasha's data analysis, HPCA+PSD95

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)
library(gridExtra)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_06_19_Sasha_HPCA-PSD')


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


df.selected <- df.full %>%
  filter(app_factor %in% c('0.5','2.5', '5', '10', '20', '30')) %>%
  mutate(rel_time = index - 30)


spines_id <- df.selected %>%
  filter(lab_id != 'shaft') %>%
  droplevels() %>%
  mutate(roi_id = interaction(roi, id, sep = '_')) %>%
  group_by(roi_id) %>%
  summarise(has_all = n_distinct(lab_id) == 2) %>%
  filter(has_all) %>%
  droplevels() %>%
  pull(roi_id)

df.spines <- df.selected %>%
  filter(lab_id != 'shaft') %>%
  mutate(roi_id = interaction(roi, id, sep = '_')) %>%
  filter(roi_id %in% spines_id) %>%
  droplevels()


##### EXPLORATORY ANALYSIS #####
df.summary <- df.spines %>%
  filter(dist_group != '0') %>%
  group_by(lab_id, app_factor, dist_group) %>%
  summarise(n_cell = n_distinct(id), n_roi = n_distinct(roi),
            max = max(df))
  
ggplot(data = df.full %>% filter(dist_group != '0'), aes(x = df, fill = dist_group)) +
  geom_density(alpha = .5) +
  facet_wrap(~interaction(lab_id, app_factor))


##### AMPLITUDE VS APP TIME #####
df.max <- df.spines %>%
  filter(dist_group != '0') %>%
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

ggplot(data = df.max, aes(x = app, y = max_df, color = lab_id, group = lab_id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.75) +
  stat_summary(fun = median,
               geom = 'point', size = 2) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'errorbar', linewidth = 0.5, width=0.05) +
  scale_x_continuous(trans = "log10", breaks = c(0.5, 2.5, 5, 10, 20, 60, 120)) +
  scale_y_continuous(breaks = seq(0,10,0.1), limits = c(-0.1,0.5)) +
  facet_wrap(~dist_group, nrow = 3) +
  theme_minimal_hgrid()



##### AMPLITUDE 20 ANALYSIS #####
df.20 <-  df.full %>%
  filter(app_factor == '0.5') %>%
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
              linewidth = 0, alpha = .3) +
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
  

##### LM AGREGATE PLOT #####
# selected frames lm
time.points <- c(3, 15, 18, 20)

ggplot(data = df.spines %>%
              filter(rel_time %in% time.points, lab_id != 'shaft', app_factor == '20') %>%
              mutate(lab_dist = interaction(lab_id, dist_group, sep = '_')),
       aes(x = dist_um, y = df,
           color = lab_id, shape = lab_id, fill = lab_id, group = lab_dist)) +
  geom_vline(xintercept = c(25, 50), linetype = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = .1) +
  # geom_line(aes(group = roi_id), alpha = .1) +
  geom_smooth(method = 'lm') +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  scale_color_manual(values = c('oreol' = 'coral', 'psd' = 'magenta3')) +
  scale_fill_manual(values = c('oreol' = 'coral', 'psd' = 'magenta3')) +
  facet_wrap(~as.factor(rel_time), ncol = length(time.points)) +  # ~interaction(as.factor(rel_time), app_factor, sep = '_')
  labs(title = 'Linear fit for selected frames (app. 20 s)',
       x = 'Distance (μm)',
       y = 'ΔF/F0',
       color = 'Region',
       fill = 'Region',
       shape = 'Region') +
  theme_minimal()



##### ROI STAT #####

# df.by.roi <- df.selected %>%
#   filter(lab_id != 'shaft', dist_group %in% c('II', 'III')) %>%
#   mutate(roi_id = interaction(roi, id, sep = '_')) %>%
#   select(roi_id, df, lab_id, rel_time, app_factor, dist_group) %>%
#   pivot_wider(names_from = lab_id, values_from = df) %>%
#   drop_na() %>%
#   group_by(roi_id) %>%
#   mutate(filter_group = ifelse(((rel_time == 2) & (oreol > psd)), TRUE, FALSE),
#          rise_group = ifelse(!all(filter_group == FALSE), 'up', 'down')) %>%
#   select(-filter_group) %>%
#   droplevels() %>%
#   ungroup()

df.by.roi <- df.spines %>%
  filter(dist_group %in% c('II', 'III')) %>%
  select(roi_id, df, lab_id, rel_time, app_factor, dist_group) %>%
  pivot_wider(names_from = lab_id, values_from = df) %>%
  drop_na() %>%
  group_by(roi_id) %>%
  mutate(filter_group = ifelse(((rel_time == 2) & (oreol > psd)), TRUE, FALSE),
         rise_group = ifelse(!all(filter_group == FALSE), 'up', 'down')) %>%
  select(-filter_group) %>%
  droplevels() %>%
  ungroup()


track.time <- seq(0,30)
track.end <- seq(19,35)

df.avg.roi <- df.by.roi %>%
  filter(app_factor %in% c('20'), dist_group != 'I') %>%  # 
  group_by(rel_time) %>%
  summarise(psd_mean = mean(psd, na.rm = TRUE),
            oreol_mean = mean(oreol, na.rm = TRUE),
            psd_se = sd(psd, na.rm = TRUE) / sqrt(n()),
            oreol_se = sd(oreol, na.rm = TRUE) / sqrt(n()),
            n = n(),
            .groups = 'drop')

# ggplot(df.by.roi %>% filter(rel_time %in% track.time),
#        aes(x = psd, y = oreol, color = rel_time)) +

#   geom_point(aes(group = roi_id), alpha = .15) +

# ggplot() +
# geom_point(data = df.by.roi %>% filter(rel_time %in% track.time),
#            aes(x = psd, y = oreol, color = rel_time), alpha = .15)


ggplot() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  # geom_line(data = df.avg.roi %>% filter(rel_time %in% track.time),  # up
  #           aes(x = psd_mean, y = oreol_mean, color = rel_time)) +
  # geom_segment(data = df.avg.roi %>% filter(rel_time %in% track.time),
  #              aes(x = psd_mean, y = oreol_mean,
  #                  xend = lead(psd_mean), yend = lead(oreol_mean),
  #                  color = rel_time),
  #              arrow = arrow(length = unit(0.3, "cm")), size = 0.75) +
  geom_path(data = df.avg.roi %>% filter(rel_time %in% track.time),
            aes(x = psd_mean, y = oreol_mean, color = rel_time)) +
  geom_point(data = df.avg.roi %>% filter(rel_time %in% track.time),
             aes(x = psd_mean, y = oreol_mean, color = rel_time)) +
  geom_linerange(data = df.avg.roi %>% filter(rel_time %in% track.time),
                aes(x = psd_mean, y = oreol_mean,
                    xmin = psd_mean - psd_se,
                    xmax = psd_mean + psd_se,
                    color = rel_time),
                size = 0.3) +
  geom_linerange(data = df.avg.roi %>% filter(rel_time %in% track.time),
                aes(x = psd_mean, y = oreol_mean,
                    ymin = oreol_mean - oreol_se,
                    ymax = oreol_mean + oreol_se,
                    color = rel_time),
                size = 0.3) +
  # geom_line(data = df.avg.roi %>% filter(rel_time %in% track.end),  # down
  #           aes(x = psd_mean, y = oreol_mean, color = rel_time)) +
  # geom_segment(data = df.avg.roi %>% filter(rel_time %in% track.end),
  #              aes(x = psd_mean, y = oreol_mean,
  #                  xend = lead(psd_mean), yend = lead(oreol_mean),
  #                  color = rel_time),
  #              arrow = arrow(length = unit(0.3, "cm")), size = 0.75) +
  # geom_point(data = df.avg.roi %>% filter(rel_time %in% track.end),
  #            aes(x = psd_mean, y = oreol_mean, color = rel_time)) +
  # geom_errorbar(data = df.avg.roi %>% filter(rel_time %in% track.end),
  #               aes(x = psd_mean, y = oreol_mean,
  #                   xmin = psd_mean - psd_se,
  #                   xmax = psd_mean + psd_se,
  #                   color = rel_time),
  #               width = 0, size = 0.3) +
  # geom_errorbar(data = df.avg.roi %>% filter(rel_time %in% track.end),
  #               aes(x = psd_mean, y = oreol_mean,
  #                   ymin = oreol_mean - oreol_se,
  #                   ymax = oreol_mean + oreol_se,
  #                   color = rel_time),
  #               width = 0, size = 0.3) +
  # scale_x_continuous(limits = c(-0.05,0.0125), breaks = seq(-1,1,0.01)) +
  # scale_y_continuous(limits = c(-0.025,0.1), breaks = seq(-1,1,0.025)) +
  scale_color_gradient2(mid = 'blue', high = "red") +
  theme_minimal()
  


##### ROI REVIEW WITH BOX #####
ggplot(data = df.spines %>% filter(app_factor == '20', rel_time == 3),
       aes(x = lab_id, y = df,
           fill = lab_id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot() +
  geom_point(alpha = .15) +
  geom_line(aes(group = roi_id), linewidth = .1) +
  facet_wrap(~dist_group, ncol = 3)

ggplot(data = df.spines %>% filter(app_factor == '20', dist_group != 'I'),
       aes(x = rel_time, y = df)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_vline(xintercept = 20, linetype = 3) +
  stat_summary(aes(color = lab_id, group = interaction(id,lab_id)),
               fun = median,
               geom = 'line', linewidth = 0.15) +
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
  theme_minimal()
  


df.by.roi.tidy <- df.by.roi %>%
  pivot_longer(cols = c('oreol', 'psd'),
               names_to = 'lab_id',
               values_to = 'df')

ggplot(data = df.by.roi.tidy %>% filter(app_factor == '20',
                                        rel_time == 7),
       aes(x = lab_id, y = df, fill = lab_id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(alpha = .5) +
  geom_point(aes(group = roi_id), alpha = .15) +
  geom_line(aes(group = roi_id), linewidth = .1)



##### 20 LINE STAT #####
df.for.stat <- df.20.2_5 %>%
  mutate(rel_time = index - 30) %>%
  filter(rel_time %in% seq(-2,25)) %>%
  select(df, lab_id, dist_group, id, rel_time, app_factor)

df.stat <- df.for.stat %>%
  filter(lab_id != 'psd') %>%
  group_by(dist_group, rel_time, app_factor) %>%
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
  mutate(idx = rel_time-19,
         df = if_else(df < 0, 0.000000000000001, df)) %>%
  filter(rel_time %in% tail.20,
         lab_id != 'psd') %>%  # select(df, lab_id, roi, id, rel_time) %>%
  distinct() %>%
  droplevels() %>%
  select(roi, lab_id, id, df, idx, dist_group) %>%
  group_by(roi, lab_id, id) %>%
  mutate(filter_group = ifelse(((idx == 1) & (df > 0.2)), TRUE, FALSE)) %>%
  filter(!all(filter_group == FALSE)) %>%
  select(-filter_group) %>%
  ungroup() %>%
  nest(.by = c(roi, lab_id, id, dist_group)) %>%
  mutate(model = map(data, ~ tidy(lm(log(df) ~ idx, data = .x)))) %>%
  unnest(model) %>%
  ungroup() %>%
  select(-std.error, -statistic, -p.value) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  mutate(decay_tau = -1/idx) %>%
  filter(decay_tau > 1, decay_tau < 100)


df.lm.20 %>%
  group_by(lab_id) %>%
  kruskal_test(decay_tau ~ dist_group)

df.lm.20 %>%
  group_by(dist_group) %>%
  wilcox_test(decay_tau ~ lab_id)


ggplot(data = df.lm.20,
       aes(fill = dist_group, y = decay_tau, x = dist_group)) +
  geom_boxplot(alpha = .5) +
  facet_wrap(~lab_id)
