# NMDA ionophoresis, FRET between AP2B1-EYFP and HPCA(WT)-ECFP
# Copyright © 2024 Borys Olifirov

require(stringr)
require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2024_08_29_AP2_FRET_60|60|180_stat')

##### DATA PREPROCESSING #####
# UP MASK
df.ch0_df.up_mask <- bind_rows(read.csv('./24_05_16_09/24_05_16_09_ch0_24_05_16_09_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_05_22_06/24_05_22_06_ch0_24_05_22_06_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_01/24_06_5_01_ch0_24_06_5_01_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_02/24_06_5_02_ch0_24_06_5_01_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_04/24_06_5_04_ch0_24_06_5_04_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_06/24_06_5_06_ch0_24_06_5_06_ch0_red-green_up-labels_ΔF.csv')) %>%
                      mutate(id = as.factor(str_remove(id, '_xform_ch0')), channel = as.factor('ch0'))

df.ch3_df.up_mask <- bind_rows(read.csv('./24_05_16_09/24_05_16_09_ch3_24_05_16_09_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_05_22_06/24_05_22_06_ch3_24_05_22_06_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_01/24_06_5_01_ch3_24_06_5_01_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_02/24_06_5_02_ch3_24_06_5_01_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_04/24_06_5_04_ch3_24_06_5_04_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_06/24_06_5_06_ch3_24_06_5_06_ch0_red-green_up-labels_ΔF.csv')) %>%
                      mutate(id = as.factor(str_remove(id, '_xform_ch3')), channel = as.factor('ch3'))

df.eapp_abs.up_mask <- bind_rows(read.csv('./24_05_16_09/24_05_16_09_Eapp_24_05_16_09_ch0_red-green_up-labels_abs.csv'),
                               read.csv('./24_05_22_06/24_05_22_06_Eapp_24_05_22_06_ch0_red-green_up-labels_abs.csv'),
                               read.csv('./24_06_5_01/24_06_5_01_Eapp_24_06_5_01_ch0_red-green_up-labels_abs.csv'),
                               read.csv('./24_06_5_02/24_06_5_02_Eapp_24_06_5_01_ch0_red-green_up-labels_abs.csv'),
                               read.csv('./24_06_5_04/24_06_5_04_Eapp_24_06_5_04_ch0_red-green_up-labels_abs.csv'),
                               read.csv('./24_06_5_06/24_06_5_06_Eapp_24_06_5_06_ch0_red-green_up-labels_abs.csv')) %>%
  mutate(id = as.factor(str_remove(id, '_xform_Eapp')), channel = as.factor('Eapp'))

df.up_mask <- bind_rows(df.ch0_df.up_mask, df.ch3_df.up_mask, df.eapp_abs.up_mask) %>%
  mutate(mask = as.factor('up'))
remove(df.ch0_df.up_mask, df.ch3_df.up_mask, df.eapp_abs.up_mask)
  
# FRET MASK
df.ch0_df.fret_mask <- bind_rows(read.csv('./24_05_16_09/24_05_16_09_ch0_24_05_16_09_FRET_up-labels_ΔF.csv'),
                                 read.csv('./24_05_22_06/24_05_22_06_ch0_24_05_22_06_FRET_spine-labels_ΔF.csv'),
                                 read.csv('./24_06_5_01/24_06_5_01_ch0_24_06_5_01_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_06_5_02/24_06_5_02_ch0_24_06_5_01_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_06_5_04/24_06_5_04_ch0_24_06_5_04_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_06_5_06/24_06_5_06_ch0_24_06_5_06_Fc_norm_red-green_up-labels_ΔF.csv')) %>%
  mutate(id = as.factor(str_remove(id, '_xform_ch0')), channel = as.factor('ch0'))

df.ch3_df.fret_mask <- bind_rows(read.csv('./24_05_16_09/24_05_16_09_ch3_24_05_16_09_FRET_up-labels_ΔF.csv'),
                                 read.csv('./24_05_22_06/24_05_22_06_ch3_24_05_22_06_FRET_spine-labels_ΔF.csv'),
                                 read.csv('./24_06_5_01/24_06_5_01_ch3_24_06_5_01_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_06_5_02/24_06_5_02_ch3_24_06_5_01_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_06_5_04/24_06_5_04_ch3_24_06_5_04_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_06_5_06/24_06_5_06_ch3_24_06_5_06_Fc_norm_red-green_up-labels_ΔF.csv')) %>%
  mutate(id = as.factor(str_remove(id, '_xform_ch3')), channel = as.factor('ch3'))

df.eapp_abs.fret_mask <- bind_rows(read.csv('./24_05_16_09/24_05_16_09_Eapp_24_05_16_09_FRET_up-labels_abs.csv'),
                                 read.csv('./24_05_22_06/24_05_22_06_Eapp_24_05_22_06_FRET_spine-labels_abs.csv'),
                                 read.csv('./24_06_5_01/24_06_5_01_Eapp_24_06_5_01_Fc_norm_red-green_up-labels_abs.csv'),
                                 read.csv('./24_06_5_02/24_06_5_02_Eapp_24_06_5_01_Fc_norm_red-green_up-labels_abs.csv'),
                                 read.csv('./24_06_5_04/24_06_5_04_Eapp_24_06_5_04_Fc_norm_red-green_up-labels_abs.csv'),
                                 read.csv('./24_06_5_06/24_06_5_06_Eapp_24_06_5_06_Fc_norm_red-green_up-labels_abs.csv')) %>%
  mutate(id = as.factor(str_remove(id, '_xform_Eapp')), channel = as.factor('Eapp'))

df.fret_mask <- bind_rows(df.ch0_df.fret_mask, df.ch3_df.fret_mask, df.eapp_abs.fret_mask) %>%
  mutate(mask = as.factor('fret'))
remove(df.ch0_df.fret_mask, df.ch3_df.fret_mask, df.eapp_abs.fret_mask)

df.mask <- bind_rows(df.up_mask, df.fret_mask)
remove(df.up_mask, df.fret_mask)

df.mask <- df.mask %>%
  mutate(dist_group = cut(dist,
                          breaks = c(0, 150, 250, 10000),
                          labels = c('max', 'mid', 'min'))) %>%
  filter(index <= 29, id != '24_06_5_04')  # , , id != '24_05_16_09'


##### CTRL PLOTS #####
df.to.plot <- df.mask %>%
  filter(dist_group == 'max')

ggplot(data = df.to.plot,
       aes(x = index, y = int, color = id, group = id)) +
  stat_summary(fun = median,
               geom = 'line', size = .5) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(data = df.to.plot,
               aes(x = index, y = int, group = channel),
               color = 'black',
               fun = median,
               geom = 'point', size = 1) +
  stat_summary(data = df.to.plot,
               aes(x = index, y = int, group = channel),
               color = 'black',
               fun = median,
               geom = 'line', size = 0.4) +
  stat_summary(data = df.to.plot,
               aes(x = index, y = int, group = channel),
               color = 'black',
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', size = 0.25, width = 1) +
  annotate('rect', xmin = 6, xmax = 12, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = 'red') +
  facet_wrap(facets = vars(channel, mask), nrow = nlevels(df.mask$channel), ncol = nlevels(df.mask$mask),
             strip.position = 'right', scales = "free_y")


##### STAT PLOTS #####
df.to.stat <- df.mask

ggplot(data = df.to.stat,
       aes(x = time, y = int, color = mask, group = mask)) +
  stat_summary(fun = median,
               geom = 'line', size = .5) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', size = .5, width = .3) +
  annotate('rect', xmin = 60, xmax = 120, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = 'red') +
  facet_wrap(facets = vars(channel, dist_group), nrow = nlevels(df.to.stat$channel),ncol = nlevels(df.to.stat$dist_group),
             strip.position = 'right')  # , scales = "free_y"


##### TIMEPOINTS BOX #####
df.to.box <- df.mask %>%
  filter(index %in% c(2, 9, 13, 25),
         channel == 'Eapp',
         dist_group == 'max') %>%
  droplevels() %>%
  group_by(id, channel, index, mask) %>%
  mutate(int_med = median(int), int_iqr = IQR(int), idx = as.factor(index)) %>%
  ungroup() %>%
  select(-X, -time, -roi, -int, -dist, -dist_group, -channel, -index) %>%
  distinct()

# I vs T
df.to.box.stat.it <- df.to.box %>%
  select(-id, -int_iqr) %>%
  group_by(mask) %>%
  pairwise_wilcox_test(int_med ~ idx, p.adjust.method = "BH") %>%
  add_xy_position(fun = "mean_sd")  

ggplot(data = df.to.box,  #  %>% filter(mask == 'fret')
       aes(x = idx, y = int_med)) +
  geom_boxplot(aes(fill = idx)) +
  stat_pvalue_manual(df.to.box.stat.it, label = 'p.adj.signif',
                     hide.ns = TRUE) +
  facet_wrap(facets = vars(mask), ncol = 2) 

# I vs M
df.to.box.stat.im <- df.to.box %>%
  select(-id, -int_iqr) %>%
  group_by(idx) %>%
  wilcox_test(int_med ~ mask) %>%
  add_significance() %>%
  add_xy_position(fun = "mean_sd")  

ggplot(data = df.to.box,
       aes(x = mask, y = int_med)) +
  geom_boxplot(aes(fill = mask)) +
  stat_pvalue_manual(df.to.box.stat.im, label = 'p.signif',
                     hide.ns = TRUE) +
  facet_wrap(facets = vars(idx), ncol = 4)

