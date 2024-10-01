# NMDA ionophoresis, FRET between AP2B1-EYFP and HPCA(WT)-ECFP
# Copyright © 2024 Borys Olifirov

require(stringr)

require(dplyr)
require(tidyr)
require(purrr)
require(rstatix)

require(mixtools)
require(Rbeast)

require(ggplot2)
require(ggforce)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2024_08_29_AP2_FRET_60|60|180_stat')

font.size <- 20
font.fam <- 'Arial'
box.alpha <- 0.6

selected.mask <- 'fret'
base.indexes <- seq(0,5)
mid.indexes.f <- seq(10,12)
mid.indexes.h <- seq(10,12)  # seq(7,9)
end.indexes <- seq(24,28)

##### DATA PREPROCESSING #####
# UP MASK
df.ch0_df.up_mask <- bind_rows(read.csv('./24_05_16_04/24_05_16_04_ch0_24_05_16_04_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_05_16_08/24_05_16_08_ch0_24_05_16_08_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_05_16_09/24_05_16_09_ch0_24_05_16_09_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_05_22_04/24_05_22_04_ch0_24_05_22_04_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_05_22_06/24_05_22_06_ch0_24_05_22_06_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_05_22_07/24_05_22_07_ch0_24_05_22_07_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_01/24_06_5_01_ch0_24_06_5_01_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_02/24_06_5_02_ch0_24_06_5_01_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_04/24_06_5_04_ch0_24_06_5_04_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_06/24_06_5_06_ch0_24_06_5_06_ch0_red-green_up-labels_ΔF.csv')) %>%
                      mutate(id = as.factor(str_remove(id, '_xform_ch0')),
                             channel = as.factor('ch0'),
                             mask = as.factor('up'),
                             int_val = as.factor('df')) %>%
                      select(-X)
df.ch0_abs.up_mask <- bind_rows(read.csv('./24_05_16_04/24_05_16_04_ch0_24_05_16_04_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_05_16_08/24_05_16_08_ch0_24_05_16_08_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_05_16_09/24_05_16_09_ch0_24_05_16_09_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_05_22_04/24_05_22_04_ch0_24_05_22_04_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_05_22_06/24_05_22_06_ch0_24_05_22_06_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_05_22_07/24_05_22_07_ch0_24_05_22_07_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_06_5_01/24_06_5_01_ch0_24_06_5_01_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_06_5_02/24_06_5_02_ch0_24_06_5_01_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_06_5_04/24_06_5_04_ch0_24_06_5_04_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_06_5_06/24_06_5_06_ch0_24_06_5_06_ch0_red-green_up-labels_abs.csv')) %>%
                     mutate(id = as.factor(str_remove(id, '_xform_ch0')),
                            channel = as.factor('ch0'),
                            mask = as.factor('up'),
                            int_val = as.factor('abs')) %>%
                     select(-X)

df.ch3_df.up_mask <- bind_rows(read.csv('./24_05_16_04/24_05_16_04_ch3_24_05_16_04_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_05_16_08/24_05_16_08_ch3_24_05_16_08_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_05_16_09/24_05_16_09_ch3_24_05_16_09_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_05_22_04/24_05_22_04_ch3_24_05_22_04_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_05_22_06/24_05_22_06_ch3_24_05_22_06_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_05_22_07/24_05_22_07_ch3_24_05_22_07_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_01/24_06_5_01_ch3_24_06_5_01_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_02/24_06_5_02_ch3_24_06_5_01_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_04/24_06_5_04_ch3_24_06_5_04_ch0_red-green_up-labels_ΔF.csv'),
                               read.csv('./24_06_5_06/24_06_5_06_ch3_24_06_5_06_ch0_red-green_up-labels_ΔF.csv')) %>%
                     mutate(id = as.factor(str_remove(id, '_xform_ch3')),
                           channel = as.factor('ch3'),
                           mask = as.factor('up'),
                           int_val = as.factor('df')) %>%
                     select(-X)
df.ch3_abs.up_mask <- bind_rows(read.csv('./24_05_16_04/24_05_16_04_ch3_24_05_16_04_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_05_16_08/24_05_16_08_ch3_24_05_16_08_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_05_16_09/24_05_16_09_ch3_24_05_16_09_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_05_22_04/24_05_22_04_ch3_24_05_22_04_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_05_22_06/24_05_22_06_ch3_24_05_22_06_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_05_22_07/24_05_22_07_ch3_24_05_22_07_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_06_5_01/24_06_5_01_ch3_24_06_5_01_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_06_5_02/24_06_5_02_ch3_24_06_5_01_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_06_5_04/24_06_5_04_ch3_24_06_5_04_ch0_red-green_up-labels_abs.csv'),
                                read.csv('./24_06_5_06/24_06_5_06_ch3_24_06_5_06_ch0_red-green_up-labels_abs.csv')) %>%
                        mutate(id = as.factor(str_remove(id, '_xform_ch3')),
                               channel = as.factor('ch3'),
                               mask = as.factor('up'),
                               int_val = as.factor('abs')) %>%
                        select(-X)

df.eapp_abs.up_mask <- bind_rows(read.csv('./24_05_16_04/24_05_16_04_Eapp_24_05_16_04_ch0_red-green_up-labels_abs.csv'),
                                 read.csv('./24_05_16_08/24_05_16_08_Eapp_24_05_16_08_ch0_red-green_up-labels_abs.csv'),
                                 read.csv('./24_05_16_09/24_05_16_09_Eapp_24_05_16_09_ch0_red-green_up-labels_abs.csv'),
                                 read.csv('./24_05_22_04/24_05_22_04_Eapp_24_05_22_04_ch0_red-green_up-labels_abs.csv'),
                                 read.csv('./24_05_22_06/24_05_22_06_Eapp_24_05_22_06_ch0_red-green_up-labels_abs.csv'),
                                 read.csv('./24_05_22_07/24_05_22_07_Eapp_24_05_22_07_ch0_red-green_up-labels_abs.csv'),
                                 read.csv('./24_06_5_01/24_06_5_01_Eapp_24_06_5_01_ch0_red-green_up-labels_abs.csv'),
                                 read.csv('./24_06_5_02/24_06_5_02_Eapp_24_06_5_01_ch0_red-green_up-labels_abs.csv'),
                                 read.csv('./24_06_5_04/24_06_5_04_Eapp_24_06_5_04_ch0_red-green_up-labels_abs.csv'),
                                 read.csv('./24_06_5_06/24_06_5_06_Eapp_24_06_5_06_ch0_red-green_up-labels_abs.csv')) %>%
                       mutate(id = as.factor(str_remove(id, '_xform_Eapp')),
                              channel = as.factor('Eapp'),
                              mask = as.factor('up'),
                              int_val = as.factor('abs')) %>%
                       select(-X)

df.up_mask <- bind_rows(df.ch0_df.up_mask,
                        df.ch0_abs.up_mask,
                        df.ch3_df.up_mask,
                        df.ch3_abs.up_mask,
                        df.eapp_abs.up_mask)
remove(df.ch0_df.up_mask,
       df.ch0_abs.up_mask,
       df.ch3_df.up_mask,
       df.ch3_abs.up_mask,
       df.eapp_abs.up_mask)
  
# FRET MASK
df.ch0_df.fret_mask <- bind_rows(read.csv('./24_05_16_04/24_05_16_04_ch0_24_05_16_04_FRET_up-labels_ΔF.csv'),
                                 read.csv('./24_05_16_08/24_05_16_08_ch0_24_05_16_08_FRET_up-labels_ΔF.csv'),
                                 read.csv('./24_05_16_09/24_05_16_09_ch0_24_05_16_09_FRET_up-labels_ΔF.csv'),
                                 read.csv('./24_05_22_04/24_05_22_04_ch0_24_05_22_04_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_05_22_06/24_05_22_06_ch0_24_05_22_06_FRET_spine-labels_ΔF.csv'),
                                 read.csv('./24_05_22_07/24_05_22_07_ch0_24_05_22_07_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_06_5_01/24_06_5_01_ch0_24_06_5_01_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_06_5_02/24_06_5_02_ch0_24_06_5_01_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_06_5_04/24_06_5_04_ch0_24_06_5_04_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_06_5_06/24_06_5_06_ch0_24_06_5_06_Fc_norm_red-green_up-labels_ΔF.csv')) %>%
                      mutate(id = as.factor(str_remove(id, '_xform_ch0')),
                             channel = as.factor('ch0'),
                             mask = as.factor('fret'),
                             int_val = as.factor('df')) %>%
                      select(-X)
df.ch0_abs.fret_mask <- bind_rows(read.csv('./24_05_16_04/24_05_16_04_ch0_24_05_16_04_FRET_up-labels_abs.csv'),
                                  read.csv('./24_05_16_08/24_05_16_08_ch0_24_05_16_08_FRET_up-labels_abs.csv'),
                                  read.csv('./24_05_16_09/24_05_16_09_ch0_24_05_16_09_FRET_up-labels_abs.csv'),
                                  read.csv('./24_05_22_04/24_05_22_04_ch0_24_05_22_04_Fc_norm_red-green_up-labels_abs.csv'),
                                  read.csv('./24_05_22_06/24_05_22_06_ch0_24_05_22_06_FRET_spine-labels_abs.csv'),
                                  read.csv('./24_05_22_07/24_05_22_07_ch0_24_05_22_07_Fc_norm_red-green_up-labels_abs.csv'),
                                  read.csv('./24_06_5_01/24_06_5_01_ch0_24_06_5_01_Fc_norm_red-green_up-labels_abs.csv'),
                                  read.csv('./24_06_5_02/24_06_5_02_ch0_24_06_5_01_Fc_norm_red-green_up-labels_abs.csv'),
                                  read.csv('./24_06_5_04/24_06_5_04_ch0_24_06_5_04_Fc_norm_red-green_up-labels_abs.csv'),
                                  read.csv('./24_06_5_06/24_06_5_06_ch0_24_06_5_06_Fc_norm_red-green_up-labels_abs.csv')) %>%
                        mutate(id = as.factor(str_remove(id, '_xform_ch0')),
                               channel = as.factor('ch0'),
                               mask = as.factor('fret'),
                               int_val = as.factor('abs')) %>%
                        select(-X)

df.ch3_df.fret_mask <- bind_rows(read.csv('./24_05_16_04/24_05_16_04_ch3_24_05_16_04_FRET_up-labels_ΔF.csv'),
                                 read.csv('./24_05_16_08/24_05_16_08_ch3_24_05_16_08_FRET_up-labels_ΔF.csv'),
                                 read.csv('./24_05_16_09/24_05_16_09_ch3_24_05_16_09_FRET_up-labels_ΔF.csv'),
                                 read.csv('./24_05_22_04/24_05_22_04_ch3_24_05_22_04_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_05_22_06/24_05_22_06_ch3_24_05_22_06_FRET_spine-labels_ΔF.csv'),
                                 read.csv('./24_05_22_07/24_05_22_07_ch3_24_05_22_07_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_06_5_01/24_06_5_01_ch3_24_06_5_01_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_06_5_02/24_06_5_02_ch3_24_06_5_01_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_06_5_04/24_06_5_04_ch3_24_06_5_04_Fc_norm_red-green_up-labels_ΔF.csv'),
                                 read.csv('./24_06_5_06/24_06_5_06_ch3_24_06_5_06_Fc_norm_red-green_up-labels_ΔF.csv')) %>%
                        mutate(id = as.factor(str_remove(id, '_xform_ch3')),
                               channel = as.factor('ch3'),
                               mask = as.factor('fret'),
                               int_val = as.factor('df')) %>%
                        select(-X)
df.ch3_abs.fret_mask <- bind_rows(read.csv('./24_05_16_04/24_05_16_04_ch3_24_05_16_04_FRET_up-labels_abs.csv'),
                                  read.csv('./24_05_16_08/24_05_16_08_ch3_24_05_16_08_FRET_up-labels_abs.csv'),
                                  read.csv('./24_05_16_09/24_05_16_09_ch3_24_05_16_09_FRET_up-labels_abs.csv'),
                                  read.csv('./24_05_22_04/24_05_22_04_ch3_24_05_22_04_Fc_norm_red-green_up-labels_abs.csv'),
                                  read.csv('./24_05_22_06/24_05_22_06_ch3_24_05_22_06_FRET_spine-labels_abs.csv'),
                                  read.csv('./24_05_22_07/24_05_22_07_ch3_24_05_22_07_Fc_norm_red-green_up-labels_abs.csv'),
                                  read.csv('./24_06_5_01/24_06_5_01_ch3_24_06_5_01_Fc_norm_red-green_up-labels_abs.csv'),
                                  read.csv('./24_06_5_02/24_06_5_02_ch3_24_06_5_01_Fc_norm_red-green_up-labels_abs.csv'),
                                  read.csv('./24_06_5_04/24_06_5_04_ch3_24_06_5_04_Fc_norm_red-green_up-labels_abs.csv'),
                                  read.csv('./24_06_5_06/24_06_5_06_ch3_24_06_5_06_Fc_norm_red-green_up-labels_abs.csv')) %>%
                        mutate(id = as.factor(str_remove(id, '_xform_ch3')),
                               channel = as.factor('ch3'),
                               mask = as.factor('fret'),
                               int_val = as.factor('abs')) %>%
                        select(-X)

df.eapp_abs.fret_mask <- bind_rows(read.csv('./24_05_16_04/24_05_16_04_Eapp_24_05_16_04_FRET_up-labels_abs.csv'),
                                   read.csv('./24_05_16_08/24_05_16_08_Eapp_24_05_16_08_FRET_up-labels_abs.csv'),
                                   read.csv('./24_05_16_09/24_05_16_09_Eapp_24_05_16_09_FRET_up-labels_abs.csv'),
                                   read.csv('./24_05_22_04/24_05_22_04_Eapp_24_05_22_04_Fc_norm_red-green_up-labels_abs.csv'),
                                   read.csv('./24_05_22_06/24_05_22_06_Eapp_24_05_22_06_FRET_spine-labels_abs.csv'),
                                   read.csv('./24_05_22_07/24_05_22_07_Eapp_24_05_22_07_Fc_norm_red-green_up-labels_abs.csv'),
                                   read.csv('./24_06_5_01/24_06_5_01_Eapp_24_06_5_01_Fc_norm_red-green_up-labels_abs.csv'),
                                   read.csv('./24_06_5_02/24_06_5_02_Eapp_24_06_5_01_Fc_norm_red-green_up-labels_abs.csv'),
                                   read.csv('./24_06_5_04/24_06_5_04_Eapp_24_06_5_04_Fc_norm_red-green_up-labels_abs.csv'),
                                   read.csv('./24_06_5_06/24_06_5_06_Eapp_24_06_5_06_Fc_norm_red-green_up-labels_abs.csv')) %>%
                          mutate(id = as.factor(str_remove(id, '_xform_Eapp')),
                                 channel = as.factor('Eapp'),
                                 mask = as.factor('fret'),
                                 int_val = as.factor('abs')) %>%
                          select(-X)

df.eapp_df.fret_mask <- bind_rows(read.csv('./24_05_16_04/24_05_16_04_Eapp_24_05_16_04_FRET_up-labels_ΔF.csv'),
                                   read.csv('./24_05_16_08/24_05_16_08_Eapp_24_05_16_08_FRET_up-labels_ΔF.csv'),
                                   read.csv('./24_05_16_09/24_05_16_09_Eapp_24_05_16_09_FRET_up-labels_ΔF.csv'),
                                   read.csv('./24_05_22_04/24_05_22_04_Eapp_24_05_22_04_Fc_norm_red-green_up-labels_ΔF.csv'),
                                   read.csv('./24_05_22_06/24_05_22_06_Eapp_24_05_22_06_FRET_spine-labels_abs.csv'),
                                   read.csv('./24_05_22_07/24_05_22_07_Eapp_24_05_22_07_Fc_norm_red-green_up-labels_ΔF.csv'),
                                   read.csv('./24_06_5_01/24_06_5_01_Eapp_24_06_5_01_Fc_norm_red-green_up-labels_ΔF.csv'),
                                   read.csv('./24_06_5_02/24_06_5_02_Eapp_24_06_5_01_Fc_norm_red-green_up-labels_ΔF.csv'),
                                   read.csv('./24_06_5_04/24_06_5_04_Eapp_24_06_5_04_Fc_norm_red-green_up-labels_ΔF.csv'),
                                   read.csv('./24_06_5_06/24_06_5_06_Eapp_24_06_5_06_Fc_norm_red-green_up-labels_ΔF.csv')) %>%
  mutate(id = as.factor(str_remove(id, '_xform_Eapp')),
         channel = as.factor('Eapp'),
         mask = as.factor('fret'),
         int_val = as.factor('df')) %>%
  select(-X)

df.fret_mask <- bind_rows(df.ch0_df.fret_mask,
                          df.ch0_abs.fret_mask,
                          df.ch3_df.fret_mask,
                          df.ch3_abs.fret_mask,
                          df.eapp_abs.fret_mask)
remove(df.ch0_df.fret_mask,
       df.ch0_abs.fret_mask,
       df.ch3_df.fret_mask,
       df.ch3_abs.fret_mask,
       df.eapp_abs.fret_mask)

df.mask <- bind_rows(df.up_mask, df.fret_mask)
remove(df.up_mask, df.fret_mask)

df.mask <- df.mask %>%
  mutate(dist_group = cut(dist,
                          breaks = c(0, 150, 250, 10000),
                          labels = c('max', 'mid', 'min'))) %>%
  filter(index <= 29, id != '24_05_16_08')  # , , id != '24_05_16_09' , id != '24_06_5_04'

# write.csv(df.mask, 'df_full.csv')


##### CTRL PLOTS #####
df.to.plot <- df.mask %>%
  filter(dist_group == 'max') %>%
  filter((channel %in% c('ch0', 'ch3') & int_val == 'df') |
         (channel == 'Eapp' & int_val == 'abs')) %>%
  droplevels()

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
remove(df.to.plot)

##### LOW vs HIGH PROFILES #####


# base.indexes <- c(0,1,2,3,4,5,6)
# mid.indexes <- c(12,13,14,15)
# end.indexes <- c(26,27,28)

# selected frames
df.fret.plot <- df.mask %>%
  filter(channel == 'Eapp',
         int_val == 'abs',
         mask == selected.mask) %>%
  mutate(roi = as.factor(roi),
         cell_id = id) %>%
  droplevels() %>%
  unite('roi_id', id:roi, sep = '-') %>%
  select(roi_id, int, index, time, mask, cell_id)

ggplot(data = df.fret.plot,
       aes(x = index, y = int, color = roi_id, group = roi_id)) +
  stat_summary(fun = median,
               geom = 'line', size = .5) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(data = df.fret.plot,
               aes(x = index, y = int, group = cell_id),
               color = 'black',
               fun = median,
               geom = 'point', size = 0.75) +
  stat_summary(data = df.fret.plot,
               aes(x = index, y = int, group = cell_id),
               color = 'black',
               fun = median,
               geom = 'line', size = 0.3) +
  stat_summary(data = df.fret.plot,
               aes(x = index, y = int, group = cell_id),
               color = 'black',
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', size = 0.15, width = 0.75) +
  annotate('rect', xmin = 6, xmax = 12, ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = 'red') +
  annotate('rect',
           xmin = base.indexes[1], xmax = rev(base.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  annotate('rect',
           xmin = mid.indexes.f[1], xmax = rev(mid.indexes.f)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  annotate('rect',
           xmin = end.indexes[1], xmax = rev(end.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  labs(title = 'FRET profiles in individual ROIs, spines',
       caption = 'Red rect - NMDA app., black rect - time intervals for analysis',
       y = expression(E[app]),
       x = 'Frame idx') +
  scale_fill_manual(values = rainbow(length(df.fret.plot$roi_id))) +
  theme_classic() +
  theme(legend.position="none") +
  facet_wrap(~cell_id)

df.hpca.plot <- df.mask %>%
  filter(channel == 'ch0',
         int_val == 'df',
         mask == selected.mask) %>%
  mutate(roi = as.factor(roi),
         cell_id = id) %>%
  droplevels() %>%
  unite('roi_id', id:roi, sep = '-') %>%
  select(roi_id, int, index, time, mask, cell_id)

ggplot(data = df.hpca.plot,
       aes(x = index, y = int, color = roi_id, group = roi_id)) +
  stat_summary(fun = median,
               geom = 'line', size = .5) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(data = df.hpca.plot,
               aes(x = index, y = int, group = cell_id),
               color = 'black',
               fun = median,
               geom = 'point', size = 0.75) +
  stat_summary(data = df.hpca.plot,
               aes(x = index, y = int, group = cell_id),
               color = 'black',
               fun = median,
               geom = 'line', size = 0.3) +
  stat_summary(data = df.hpca.plot,
               aes(x = index, y = int, group = cell_id),
               color = 'black',
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', size = 0.15, width = 0.75) +
  annotate('rect', xmin = 6, xmax = 12, ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = 'red') +
  annotate('rect',
           xmin = base.indexes[1], xmax = rev(base.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  annotate('rect',
           xmin = mid.indexes.h[1], xmax = rev(mid.indexes.h)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  annotate('rect',
           xmin = end.indexes[1], xmax = rev(end.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black', size = 0) +
  labs(title = 'HPCA insertions profiles in individual ROIs, spines',
       caption = 'Red rect - NMDA app., black rect - time intervals for analysis',
       y = expression(ΔF/F[0]),
       x = 'Frame idx') +
  scale_fill_manual(values = rainbow(length(df.fret.plot$roi_id))) +
  theme_classic() +
  theme(legend.position="none") +
  facet_wrap(~cell_id)


# ggsave('fret_prof.png', p1, dpi = 300)


##### L and H GMM #####

# https://www.statisticalaid.com/independent-component-analysis-ica-using-r/
# https://www.r-bloggers.com/2011/08/fitting-mixture-distributions-with-the-r-package-mixtools/ // https://stackoverflow.com/questions/52082543/curl-package-not-available-for-several-r-packages
calc.mixmdl <- function(input_vector) {
  mixmdl <- normalmixEM(input_vector)
  input_dens <- density(input_vector)
  comp1 <- dnorm(x = input_dens$x,
                 mean = mixmdl$mu[1],
                 sd = mixmdl$sigma[1]) * mixmdl$lambda[1]
  comp2 <- dnorm(x = input_dens$x,
                 mean = mixmdl$mu[2],
                 sd = mixmdl$sigma[2]) * mixmdl$lambda[2]
  return(list(data.frame(val = input_dens$x,
                         raw = input_dens$y,
                         comp1 = comp1,
                         comp2 = comp2,
                         comp_comb = comp1+comp2),
              mixmdl))
}


# calc.mixmdl.optim <- function(input_vector) {
#   comp.vals <- seq(2,10)
#   loglik.vector <- c()
#   for (comp in comp.vals) {
#     mixmdl <- normalmixEM(input_vector, k = comp)
#     loglik.vector <- append(loglik.vector, mixmdl$loglik)
#   }
#   return(loglik.vector)
# }

# calc.mixmdl3 <- function(input_vector) {
#   mixmdl <- normalmixEM(input_vector, k = 3)
#   input_dens <- density(input_vector)
#   comp1 <- dnorm(x = input_dens$x,
#                  mean = mixmdl$mu[1],
#                  sd = mixmdl$sigma[1]) * mixmdl$lambda[1]
#   comp2 <- dnorm(x = input_dens$x,
#                  mean = mixmdl$mu[2],
#                  sd = mixmdl$sigma[2]) * mixmdl$lambda[2]
#   comp3 <- dnorm(x = input_dens$x,
#                  mean = mixmdl$mu[3],
#                  sd = mixmdl$sigma[3]) * mixmdl$lambda[3]
#   return(data.frame(val = input_dens$x,
#                     raw = input_dens$y,
#                     comp1 = comp1,
#                     comp2 = comp2,
#                     comp3 = comp3,
#                     comp_comb = comp1+comp2+comp3))
# }

# bic profiles
df.gmm.optim <- bind_rows(data.frame('BIC' = c(-562.3872730565264,
                                               -548.9248686430036,
                                               -540.2812473660222,
                                               -534.7881663801178,
                                               -522.7028024160293,
                                               -509.6392893366678,
                                               -496.5527725044958,
                                               -487.71744383341616,
                                               -474.55609636508973,
                                               -461.38142725089693),
                                     'group' = as.factor('hpca_base')),
                          data.frame('BIC' = c(-42.44625226621264,
                                               -30.687724491879592,
                                               -23.719118172169296,
                                               -11.866538799947541,
                                               -2.233073800841808,
                                               10.072656270319115,
                                               16.710397718092352,
                                               26.166101672533898,
                                               37.26514642918363,
                                               43.26662210945581),
                                     'group' = as.factor('hpca_mid')),
                          data.frame('BIC' = c(-104.31193218627054,
                                               -98.11314417077497,
                                               -86.66215271265438,
                                               -81.75972081419292,
                                               -69.16733434513318,
                                               -56.11365325504556,
                                               -45.10316725844953,
                                               -37.160264743770625,
                                               -24.505149596745184,
                                               -14.171525094894577),
                                     'group' = as.factor('hpca_end')),
                          data.frame('BIC' = c(-406.07752028766777,
                                               -430.1786832145669,
                                               -430.60262901300746,
                                               -418.81525205207737,
                                               -407.563531079979,
                                               -394.8903970851102,
                                               -399.4152886875826,
                                               -386.84349245916104,
                                               -374.80243542201333,
                                               -362.60452716617715),
                                     'group' = as.factor('fret_base')),
                          data.frame('BIC' = c(-371.79158965917827,
                                               -384.40132720733305,
                                               -386.05050342658086,
                                               -373.4270002009531,
                                               -374.5023387992349,
                                               -355.3847141473434,
                                               -352.20846465988757,
                                               -340.88784781535423,
                                               -328.67362615634954,
                                               -316.0765789697767),
                                     'group' = as.factor('fret_mid')),
                          data.frame('BIC' = c(-368.4586826816587,
                                               -390.8211030869018,
                                               -383.1737634187368,
                                               -376.24180469985055,
                                               -367.8500322192769,
                                               -359.02616960850384,
                                               -348.57977152252045,
                                               -334.3279633084872,
                                               -326.4165044176185,
                                               -313.94031562272346),
                                     'group' = as.factor('fret_end')))


##### HIST DATA FRAMES #####

# base
df.fret.base <- df.mask %>%
  filter(channel == 'Eapp',
         index %in% base.indexes) %>%
  select(id, roi, int, mask, channel) %>%
  droplevels() %>%
  group_by(id, roi, mask) %>%
  mutate(int = mean(int),
         roi = as.factor(roi),
         cell_id = id) %>%
  ungroup() %>%
  distinct() %>%
  unite('roi_id', id:roi, sep = '-') %>%
  mutate(time_interval = 'base')

df.hpca.base <- df.mask %>%
  filter(channel == 'ch0', int_val == 'df',
         index %in% base.indexes) %>%
  select(id, roi, int, mask, channel) %>%
  droplevels() %>%
  group_by(id, roi, mask) %>%
  mutate(int = mean(int),
         roi = as.factor(roi),
         cell_id = id) %>%
  ungroup() %>%
  distinct() %>%
  unite('roi_id', id:roi, sep = '-') %>%
  mutate(time_interval = 'base')

bm <- calc.mixmdl(df.fret.base$int[df.fret.base$mask == selected.mask]) 
df.base_mixmdl <- bm[[1]] %>%
  mutate(channel = 'Eapp', time_interval = 'base')
base_mixmgl <- bm[[2]]
remove(bm)

# bm.h <- calc.mixmdl(df.hpca.base$int[df.hpca.base$mask == selected.mask]) 
# df.base_mixmdl.h <- bm.h[[1]] %>%
#   mutate(channel = 'ch0', time_interval = 'base')
# base_mixmgl_h <- bm.h[[2]]
# remove(bm.h)
# 
# # mid
# df.fret.mid <- df.mask %>%
#   filter(channel == 'Eapp',
#          index %in% mid.indexes.f) %>%
#   select(id, roi, int, mask, channel) %>%
#   droplevels() %>%
#   group_by(id, roi, mask) %>%
#   mutate(int = mean(int),
#          roi = as.factor(roi),
#          cell_id = id) %>%
#   ungroup() %>%
#   distinct() %>%
#   unite('roi_id', id:roi, sep = '-') %>%
#   mutate(time_interval = 'mid')
# 
# df.hpca.mid <- df.mask %>%
#   filter(channel == 'ch0', int_val == 'df',
#          index %in% mid.indexes.h) %>%
#   select(id, roi, int, mask, channel) %>%
#   droplevels() %>%
#   group_by(id, roi, mask) %>%
#   mutate(int = mean(int),
#          roi = as.factor(roi),
#          cell_id = id) %>%
#   ungroup() %>%
#   distinct() %>%
#   unite('roi_id', id:roi, sep = '-') %>%
#   mutate(time_interval = 'mid')
# 
# mm <- calc.mixmdl(df.fret.mid$int[df.fret.mid$mask == selected.mask]) 
# df.mid_mixmdl <- mm[[1]] %>%
#   mutate(channel = 'Eapp', time_interval = 'mid')
# mid_mixmgl <- mm[[2]]
# remove(mm)
# 
# mm.h <- calc.mixmdl(df.hpca.mid$int[df.hpca.mid$mask == selected.mask]) 
# df.mid_mixmdl.h <- mm.h[[1]] %>%
#   mutate(channel = 'ch0', time_interval = 'mid')
# mid_mixmgl_h <- mm.h[[2]]
# remove(mm.h)
# 
# # end
# df.fret.end <- df.mask %>%
#   filter(channel == 'Eapp',
#          index %in% end.indexes) %>%
#   select(id, roi, int, mask, channel) %>%
#   droplevels() %>%
#   group_by(id, roi, mask) %>%
#   mutate(int = mean(int),
#          roi = as.factor(roi),
#          cell_id = id) %>%
#   ungroup() %>%
#   distinct() %>%
#   unite('roi_id', id:roi, sep = '-') %>%
#   mutate(time_interval = 'end')
# 
# df.hpca.end <- df.mask %>%
#   filter(channel == 'ch0', int_val == 'df',
#          index %in% end.indexes) %>%
#   select(id, roi, int, mask, channel) %>%
#   droplevels() %>%
#   group_by(id, roi, mask) %>%
#   mutate(int = mean(int),
#          roi = as.factor(roi),
#          cell_id = id) %>%
#   ungroup() %>%
#   distinct() %>%
#   unite('roi_id', id:roi, sep = '-') %>%
#   mutate(time_interval = 'end')
# 
# em <- calc.mixmdl(df.fret.end$int[df.fret.end$mask == selected.mask]) 
# df.end_mixmdl <- em[[1]] %>%
#   mutate(channel = 'Eapp', time_interval = 'end')
# end_mixmgl <- em[[2]]
# remove(em)
# 
# em.h <- calc.mixmdl(df.hpca.end$int[df.hpca.end$mask == selected.mask]) 
# df.end_mixmdl.h <- em.h[[1]] %>%
#   mutate(channel = 'ch0', time_interval = 'end')
# end_mixmgl_h <- em.h[[2]]
# remove(em.h)

# 
# df.hist <- bind_rows(df.hpca.base, df.hpca.mid, df.hpca.end,
#                      df.fret.base, df.fret.mid, df.fret.end) %>%
#   mutate_at(c('roi_id', 'cell_id'), as.factor)
# remove(df.hpca.base, df.hpca.mid, df.hpca.end,
#        df.fret.base, df.fret.mid, df.fret.end)
# 
# df.gmm <- bind_rows(df.base_mixmdl.h, df.mid_mixmdl.h, df.end_mixmdl.h,
#                     df.base_mixmdl, df.mid_mixmdl, df.end_mixmdl) %>%
#   mutate_at(c('channel', 'time_interval'), as.factor)
# remove(df.base_mixmdl.h, df.mid_mixmdl.h, df.end_mixmdl.h,
#        df.base_mixmdl, df.mid_mixmdl, df.end_mixmdl)
# 

##### HIST PLOT #####

plot_hist <- function(df.h, df.g, df.b, mm, scale = 3,
                      bic.title = 'Ba', lims = c(-0.025, 0.1)) {
  plot.base <- ggplot() +
    geom_histogram(data = df.h,
                   aes(x = int),
                   color = 'black',
                   fill = 'white',
                   size = 1) +
    geom_line(data = df.g,
              aes(x = val, y = comp_comb / scale), size = 1.5) +
    geom_line(data = df.g,
              aes(x = val, y = comp1 / scale),
              color = 'green3', size = 1.5, alpha = .75) +
    geom_line(data = df.g,
              aes(x = val, y = comp2 / scale),
              color = 'red', size = 1.5, alpha = .75) +
    geom_vline(xintercept = mm$mu, lty = 2, size = 0.75) +
    theme_classic() +
    theme(text=element_text(size=font.size+2, family=font.fam, face="bold"),
          legend.position="none") +
    labs(y = '# ROIs',
         x = expression(E[app])) +
    scale_x_continuous(limits = lims)
  
  bic.hpca <- ggplot() +
    geom_line(data = df.b,
              aes(x = seq(1,10), y = BIC), color = 'black', size = 1) +
    geom_point(data = df.b,
               aes(x = seq(1,10), y = BIC), color = 'black', size = 2.5) +
    scale_x_continuous(breaks = seq(1,10,2)) + 
    theme_classic() +
    theme(text=element_text(size=font.size - 3, family = font.fam, face="bold"),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) +
    labs(x = '# comp.',
         y = 'BIC')
  
  ggdraw(plot.base) +
    draw_plot(bic.hpca, x = 0.6, y = 0.5, width = 0.3, height = 0.4) +
    draw_plot_label(c(bic.title),
                    c(0.12),
                    c(1),
                    size = font.size + 5)
}



# f.hist.m <- plot_hist(df.h = df.hist %>% filter(mask == selected.mask, channel == 'Eapp', time_interval == 'mid'),
#           df.g = df.gmm %>% filter(channel == 'Eapp', time_interval == 'mid'),
#           df.b = df.gmm.optim  %>% filter(group == 'fret_mid'),
#           mm = mid_mixmgl,
#           bic.title = 'Bb')
# f.hist.e <- plot_hist(df.h = df.hist %>% filter(mask == selected.mask, channel == 'Eapp', time_interval == 'end'),
#           df.g = df.gmm %>% filter(channel == 'Eapp', time_interval == 'end'),
#           df.b = df.gmm.optim  %>% filter(group == 'fret_end'),
#           mm = end_mixmgl,
#           bic.title = 'Bc')

# f.hist <- plot_grid(f.hist.b, f.hist.m, f.hist.e, ncol = 1)
# f.hist
# # hpca hist
# h.hist.b <- plot_hist(df.h = df.hist %>% filter(mask == selected.mask, channel == 'ch0', time_interval == 'base'),
#           df.g = df.gmm %>% filter(channel == 'ch0', time_interval == 'base'),
#           df.b = df.gmm.optim  %>% filter(group == 'hpca_base'),
#           mm = base_mixmgl_h,
#           scale = 0.28,
#           lims = c(-0.5,0.5),
#           bic.title = 'Ca')
# h.hist.m <- plot_hist(df.h = df.hist %>% filter(mask == selected.mask, channel == 'ch0', time_interval == 'mid'),
#           df.g = df.gmm %>% filter(channel == 'ch0', time_interval == 'mid'),
#           df.b = df.gmm.optim  %>% filter(group == 'hpca_mid'),
#           mm = mid_mixmgl_h,scale = 0.28,
#           lims = c(-0.5,0.5),
#           bic.title = 'Cb')
# h.hist.e <- plot_hist(df.h = df.hist %>% filter(mask == selected.mask, channel == 'ch0', time_interval == 'end'),
#           df.g = df.gmm %>% filter(channel == 'ch0', time_interval == 'end'),
#           df.b = df.gmm.optim  %>% filter(group == 'hpca_end'),
#           mm = end_mixmgl_h,
#           scale = 0.28,
#           lims = c(-0.5,0.5),
#           bic.title = 'Cc')
# 
# h.hist <- plot_grid(h.hist.b, h.hist.m, h.hist.e, ncol = 1)
# h.hist

# plot_grid(f.hist, h.hist, nrow = 1)

# # base hpca hist
# plot.base.h <- ggplot() +
#   geom_histogram(data = df.hpca.base %>% filter(mask == selected.mask),
#                  aes(x = int, fill = roi_id),
#                  stackgroups = TRUE, binwidth = b.width,
#                  binpositions = "all", method = "histodot") +
#   geom_line(data = df.base_mixmdl.h,
#             aes(x = val, y = comp1 / comp.scale.factor),
#             color = 'red') +
#   geom_line(data = df.base_mixmdl.h,
#             aes(x = val, y = comp2 / comp.scale.factor),
#             color = 'blue') +
#   geom_line(data = df.base_mixmdl.h,
#             aes(x = val, y = comp_comb / comp.scale.factor),
#             lty = 2) +
#   geom_vline(xintercept = base_mixmgl_h$mu, lty = 3) +
#   scale_fill_manual(values = rainbow(length(df.hpca.base$roi_id[df.hpca.base$mask == selected.mask]))) +
#   theme_classic() +
#   theme(legend.position="none",
#         axis.text.y=element_blank(), 
#         axis.ticks.y=element_blank(),
#         axis.title.y = element_blank()) +
#   labs(title = 'Base app. average HPCA insertions in individual ROIs',
#        x = expression(ΔF/F[0])) +
#   xlim(-0.5,1)


# df.for.gmm <- data.frame('fret_base' = df.fret.base$int[df.fret.base$mask == selected.mask],
#                          'hpca_base' = df.hpca.base$int[df.hpca.base$mask == selected.mask],
#                          'fret_mid' = df.fret.mid$int[df.fret.mid$mask == selected.mask],
#                          'hpca_mid' = df.hpca.mid$int[df.hpca.mid$mask == selected.mask],
#                          'fret_end' = df.fret.end$int[df.fret.end$mask == selected.mask],
#                          'hpca_end' = df.hpca.end$int[df.hpca.end$mask == selected.mask])
# write.csv(df.for.gmm, file = 'fret_for_gmm.csv')
# 
# plot(density(df.hpca.mid$int[df.hpca.mid$mask == selected.mask]))


##### TIME INTERVALS FRET #####
sites.threshold <- 0.012
index.1 <- base.indexes #  seq(0,5)
index.2 <- mid.indexes.f #  seq(10,12) #  seq(12,16)
index.3 <- end.indexes #  seq(24,28)

df.fret.sites <- df.mask %>%
  filter(channel == 'Eapp',
         int_val == 'abs',
         mask == 'fret') %>%
  mutate(roi = as.factor(roi),
         cell_id = id) %>%
  droplevels() %>%
  unite('roi_id', id:roi, sep = '-') %>%
  group_by(roi_id) %>%
  mutate(wash_tail_int = mean(int[index %in% seq(0,6)]),
         site_type = as.factor(ifelse(wash_tail_int > sites.threshold, 'Spine', 'Shaft'))) %>%
  ungroup() %>%
  select(-dist, -dist_group, -int_val, -channel, -mask)

spine_rois <- unique(df.fret.sites$roi_id[df.fret.sites$site_type == 'Spine'])

# profiles for fret, abs
plot.fret.abs <- ggplot() +
  geom_hline(yintercept = 0, lty = 2) +
  annotate('rect', xmin = 0, xmax = 60, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = 'red') +
  annotate('rect',
           xmin = (index.1[1]*10)-60, xmax = (rev(index.1)[1]*10)-60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, color = 'black', linewidth = 0) +
  annotate('rect',
           xmin = (index.2[1]*10)-60, xmax = (rev(index.2)[1]*10)-60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, color = 'black', linewidth = 0) +
  annotate('rect',
           xmin = (index.3[1]*10)-60, xmax = (rev(index.3)[1]*10)-60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, color = 'black', linewidth = 0) +
  annotate('text', label = 'I', x = -35, y = 0.06, color = 'black', size = font.size-12) +
  annotate('text', label = 'II', x = 50, y = 0.06, color = 'black', size = font.size-12) +
  annotate('text', label = 'III', x = 200, y = 0.06, color = 'black', size = font.size-12) +
  stat_summary(data = df.fret.sites,
               aes(x = time - 60, y = int, color=site_type, group = site_type),
               fun = median,
               geom = 'line', size = 0.5) +
  stat_summary(data = df.fret.sites,
               aes(x = time - 60, y = int, color = site_type, group = site_type),
               fun = median,
               geom = 'point', size = 1) +
  stat_summary(data = df.fret.sites,
               aes(x = time - 60, y = int,
                   color = site_type, fill = site_type, group = site_type),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'ribbon', size = 0, alpha = .25) +
  scale_color_manual(values = c('Spine' = 'red', 'Shaft' = 'green4')) +
  scale_fill_manual(values = c('Spine' = 'red', 'Shaft' = 'green4')) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam, face="bold")) +
  scale_x_continuous(breaks = seq(-60, 230, 20)) +
  labs(caption = 'n = 3/9/81 (cultures/cells/ROIs)',
       x = 'Time, s',
       y = expression(E[app]))

# profiles for fret, df
df.eapp_df.prof <- df.eapp_df.fret_mask %>%
  mutate(roi = as.factor(roi),
         cell_id = id) %>%
  droplevels() %>%
  unite('roi_id', id:roi, sep = '-') %>%
  mutate(roi_id = as.factor(roi_id),
         site_type = as.factor(if_else(roi_id %in% spine_rois, 'Spine', 'Shaft'))) %>%
  select(roi_id, int, index, cell_id, site_type, time) %>%
  group_by(roi_id, site_type) %>%
  mutate(time_interval = case_when(index %in% base.indexes ~ 'I',
                                   index %in% mid.indexes.f ~ 'II',
                                   index %in% end.indexes ~ 'III',
                                   .default = 'out')) %>%
  droplevels() %>%
  ungroup() %>%
  mutate(time_interval = factor(time_interval, c('I', 'II', 'III'), ordered = TRUE))

plot.fret.df <- ggplot() + geom_hline(yintercept = 0, lty = 2) +
  annotate('rect', xmin = 0, xmax = 60, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = 'red') +
  annotate('rect',
           xmin = (index.1[1]*10)-60, xmax = (rev(index.1)[1]*10)-60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, color = 'black', linewidth = 0) +
  annotate('rect',
           xmin = (index.2[1]*10)-60, xmax = (rev(index.2)[1]*10)-60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, color = 'black', linewidth = 0) +
  annotate('rect',
           xmin = (index.3[1]*10)-60, xmax = (rev(index.3)[1]*10)-60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, color = 'black', linewidth = 0) +
  annotate('text', label = 'I', x = -35, y = 7, color = 'black', size = font.size-12) +
  annotate('text', label = 'II', x = 50, y = 7, color = 'black', size = font.size-12) +
  annotate('text', label = 'III', x = 200, y = 7, color = 'black', size = font.size-12) +
  stat_summary(data = df.eapp_df.prof,
               aes(x = time - 60, y = int, color=site_type, group = site_type),
               fun = median,
               geom = 'line', size = 0.5) +
  stat_summary(data = df.eapp_df.prof,
               aes(x = time - 60, y = int, color = site_type, group = site_type),
               fun = median,
               geom = 'point', size = 1) +
  stat_summary(data = df.eapp_df.prof,
               aes(x = time - 60, y = int,
                   color = site_type, fill = site_type, group = site_type),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'ribbon', size = 0, alpha = .25) +
  scale_color_manual(values = c('Spine' = 'red', 'Shaft' = 'green4')) +
  scale_fill_manual(values = c('Spine' = 'red', 'Shaft' = 'green4')) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam, face="bold")) +
  scale_x_continuous(breaks = seq(-60, 230, 20)) +
  labs(caption = 'n = 3/9/81 (cultures/cells/ROIs)',
       x = 'Time, s',
       y = expression(ΔE[app]/E[app[0]]))

# # plat det
# df.fret.sites.med <- df.fret.sites %>%
#   group_by(site_type, index) %>%
#   mutate(int_med = median(int)) %>%
#   select(-roi_id, -cell_id, -wash_tail_int, -int) %>%
#   distinct() %>%
#   ungroup()
# 
# o1 = beast(df.fret.sites.med$int_med[df.fret.sites.med$site_type == 'Spine'],
#            season = 'none', method='bic')
# df.fret.sites.beast.spine <- data.frame('Time' = seq(0, length(o1$trend$slpSgnPosPr)-1) * 10,
#                                         'Pos.' = o1$trend$slpSgnPosPr,
#                                         'Zero' = o1$trend$slpSgnZeroPr,
#                                         'Neg.' = rep(1, length(o1$trend$slpSgnPosPr)) - (o1$trend$slpSgnPosPr+o1$trend$slpSgnZeroPr)) %>%
#   pivot_longer(cols = 'Pos.':'Neg.', names_to = 'SlopeSign', values_to = 'p') %>%
#   mutate(SlopeSign = factor(SlopeSign, c('Pos.', 'Zero', 'Neg.'), ordered = TRUE))
# ggplot() +
#   geom_area(data = df.fret.sites.beast.spine,
#             aes(x = Time, y = p, fill = SlopeSign), alpha = .6) +
#   annotate('text', label = 'NMDA app.', x = 90, y = 1.03, color = 'red') +
#   annotate('rect', xmin = 60, xmax = 120, ymin = 0, ymax = 1.06,
#            alpha = 0.2, fill = 'red') +
#   geom_hline(yintercept = 0.82, lty = 2) +
#   scale_fill_manual(values = c('Pos.' = 'red3',
#                                'Zero' = 'yellow3',
#                                'Neg.' =  'blue3')) +
#   scale_x_continuous(limits = c(0,280), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(-0.002,1.06), expand = c(0, 0)) +
#   theme_classic() +
#   labs(title = 'Probability of slope sign for spines sites',
#        fill = 'Slope sing',
#        x = 'Time, s',
#        y = 'Probability')  
# 
# o2 = beast(df.fret.sites.med$int_med[df.fret.sites.med$site_type == 'Shaft'],
#            season = 'none', method='bic')
# df.fret.sites.beast.shaft <- data.frame('Time' = seq(0, length(o2$trend$slpSgnPosPr)-1) * 10,
#                                         'Pos.' = o2$trend$slpSgnPosPr,
#                                         'Zero' = o2$trend$slpSgnZeroPr,
#                                         'Neg.' = rep(1, length(o2$trend$slpSgnPosPr)) - (o2$trend$slpSgnPosPr+o2$trend$slpSgnZeroPr)) %>%
#   pivot_longer(cols = 'Pos.':'Neg.', names_to = 'SlopeSign', values_to = 'p') %>%
#   mutate(SlopeSign = factor(SlopeSign, c('Pos.', 'Zero', 'Neg.'), ordered = TRUE))
# ggplot() +
#   geom_area(data = df.fret.sites.beast.shaft,
#             aes(x = Time, y = p, fill = SlopeSign), alpha = .6) +
#   annotate('text', label = 'NMDA app.', x = 90, y = 1.03, color = 'red') +
#   annotate('rect', xmin = 60, xmax = 120, ymin = 0, ymax = 1.06,
#            alpha = 0.2, fill = 'red') +
#   geom_hline(yintercept = 1, lty = 2) +
#   scale_fill_manual(values = c('Pos.' = 'red3',
#                                'Zero' = 'yellow3',
#                                'Neg.' =  'blue3')) +
#   scale_x_continuous(limits = c(0,280), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(-0.002,1.06), expand = c(0, 0)) +
#   theme_classic() +
#   labs(title = 'Probability of slope sign for shaft sites',
#        fill = 'Slope sing',
#        x = 'Time, s',
#        y = 'Probability')  

# time intervals
df.fret.sites.avg <- df.fret.sites %>%
  select(-wash_tail_int, -time) %>%
  group_by(roi_id, site_type) %>%
  mutate(time_interval = case_when(index %in% base.indexes ~ 'I',
                                   index %in% mid.indexes.f ~ 'II',
                                   index %in% end.indexes ~ 'III',
                                   .default = 'out')) %>%
  filter(time_interval != 'out') %>%
  droplevels() %>%
  ungroup() %>%
  mutate(time_interval = factor(time_interval, c('I', 'II', 'III'), ordered = TRUE)) %>%
  group_by(roi_id, site_type, time_interval) %>%
  mutate(int_interval = median(int)) %>%
  select(-index, -int) %>%
  ungroup() %>%
  distinct()

# base vs zero
df.base.zero.stat <- df.fret.sites.avg %>%
  filter(time_interval == 'I') %>%
  select(-time_interval) %>%
  group_by(site_type) %>%
  wilcox_test(int_interval ~ 1, mu = 0) %>%
  add_significance() %>%
  mutate(y.position = c(0.082,0.035), group2 = c(1,1))

plot.box.zero.abs <- ggplot() + geom_boxplot(data = df.fret.sites.avg %>% filter(time_interval == 'I'),
               aes(x = time_interval,
                   y = int_interval,
                   fill = site_type),
               alpha = box.alpha) +
  geom_point(data = df.fret.sites.avg %>% filter(time_interval == 'I'),
             aes(x = time_interval,
                 y = int_interval),
             size=2, shape=21) +
  stat_pvalue_manual(df.base.zero.stat, label = 'p.signif',
                     hide.ns = TRUE, remove.bracket = TRUE, size = font.size - 12) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c('Spine' = 'red', 'Shaft' = 'green4')) +
  scale_y_continuous(limits = c(0, 0.085), breaks = seq(-1,1,0.025)) +
  facet_wrap(~site_type) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam, face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')) +
  labs(x = 'Time interval',
       y = expression(E[app]))

# time vs int
df.time.interval.stat <- df.fret.sites.avg %>%
  group_by(site_type) %>%
  pairwise_wilcox_test(int_interval ~ time_interval, p.adjust.method = 'BH') %>%
  add_significance() %>%
  add_xy_position(fun = "max") 

plot.box.ti.abs <- ggplot() +
  geom_boxplot(data = df.fret.sites.avg,
               aes(x = time_interval,
                   y = int_interval,
                   fill = site_type),
               alpha = box.alpha) +
  geom_point(data = df.fret.sites.avg,
             aes(x = time_interval,
                 y = int_interval),
             size=2, shape=21) +
  geom_line(data = df.fret.sites.avg,
            aes(x = time_interval,
                y = int_interval,
                group = roi_id),
            size = .25, lty = 2, alpha = .5) +
  stat_pvalue_manual(df.time.interval.stat, label = 'p.adj.signif',
                     hide.ns = TRUE, size = font.size - 12) +
  scale_fill_manual(values = c('Spine' = 'red', 'Shaft' = 'green4')) +
  facet_wrap(~site_type) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam, face="bold")) +
  labs(x = 'Time interval',
       y = expression(E[app]))


# site vs int
df.site.type.stat <- df.fret.sites.avg %>%
  group_by(time_interval) %>%
  wilcox_test(int_interval ~ site_type) %>%
  add_significance() %>%
  add_xy_position(fun = "max") 

plot.box.si.abs <-ggplot() +
  geom_boxplot(data = df.fret.sites.avg,
               aes(x = site_type,
                   y = int_interval,
                   fill = site_type),
               alpha = box.alpha) +
  geom_point(data = df.fret.sites.avg,
             aes(x = site_type,
                 y = int_interval),
             size=2, shape=21) +
  stat_pvalue_manual(df.site.type.stat, label = 'p.signif',
                     hide.ns = TRUE, size = font.size - 12) +
  scale_fill_manual(values = c('Spine' = 'red', 'Shaft' = 'green4')) +
  facet_wrap(~time_interval) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam, face="bold"),
        axis.text.x = element_text(angle = -90, hjust = 1)) +
  labs(x = 'Site type',
       y = expression(E[app]))



# time vs int df
df.time.interval.df <- df.eapp_df.prof %>%
  filter(time_interval != 'out') %>%
  group_by(roi_id, site_type, time_interval) %>%
  mutate(int_interval = median(int)) %>%
  ungroup() %>%
  filter(int_interval < 200) %>%
  mutate(site_type = factor(site_type, c('Spine', 'Shaft'), ordered = TRUE))

df.time.interval.df.stat <- df.time.interval.df %>%
  group_by(site_type) %>%
  pairwise_wilcox_test(int_interval ~ time_interval, p.adjust.method = 'BH') %>%
  add_significance() %>%
  add_xy_position(fun = "max") 


plot.box.ti.df <- ggplot() +
  geom_boxplot(data = df.time.interval.df,
               aes(x = time_interval,
                   y = int_interval,
                   fill = site_type),
               alpha = box.alpha) +
  geom_point(data = df.time.interval.df,
             aes(x = time_interval,
                 y = int_interval),
             size=2, shape=21) +
  geom_line(data = df.time.interval.df,
            aes(x = time_interval,
                y = int_interval,
                group = roi_id),
            size = .25, lty = 2, alpha = .5) +
  stat_pvalue_manual(df.time.interval.df.stat, label = 'p.adj.signif',
                     hide.ns = TRUE, size = font.size - 12) +
  scale_fill_manual(values = c('Spine' = 'red', 'Shaft' = 'green4')) +
  facet_wrap(~site_type, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam, face="bold")) +
  labs(x = 'Time interval',
       y = expression(ΔE[app]/E[app[0]]))


# site vs int df
df.site.type.df.stat <- df.time.interval.df %>%
  group_by(time_interval) %>%
  wilcox_test(int_interval ~ site_type) %>%
  add_significance() %>%
  add_xy_position(fun = "max") 

plot.box.si.df <- ggplot() +
  geom_boxplot(data = df.time.interval.df ,
               aes(x = site_type,
                   y = int_interval,
                   fill = site_type),
               alpha = box.alpha) +
  geom_point(data = df.time.interval.df ,
             aes(x = site_type,
                 y = int_interval),
             size=2, shape=21) +
  stat_pvalue_manual(df.site.type.df.stat, label = 'p.signif',
                     hide.ns = TRUE, size = font.size - 12) +
  scale_fill_manual(values = c('Spine' = 'red', 'Shaft' = 'green4')) +
  facet_wrap(~time_interval, scales = "free") +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam, face="bold"),
        axis.text.x = element_text(angle = -90, hjust = 1)) +
  labs(x = 'Site type',
       y = expression(ΔE[app]/E[app[0]]))


##### TIME INTERVALS HPCA #####
df.hpca.sites <- df.mask %>%
  filter(channel == 'ch0',
         int_val == 'df',
         mask == 'fret') %>%
  mutate(roi = as.factor(roi),
         cell_id = id) %>%
  droplevels() %>%
  unite('roi_id', id:roi, sep = '-') %>%
  mutate(site_type = df.fret.sites$site_type) %>%
  ungroup() %>%
  select(-dist, -dist_group, -int_val, -channel, -mask)

# profiles
plot.hpca.df <- ggplot() +
  annotate('rect', xmin = 0, xmax = 60, ymin = -Inf, ymax = Inf,
         alpha = 0.15, fill = 'red') +
  annotate('rect',
           xmin = (index.1[1]*10)-60, xmax = (rev(index.1)[1]*10)-60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, color = 'black', linewidth = 0) +
  annotate('rect',
           xmin = (index.2[1]*10)-60, xmax = (rev(index.2)[1]*10)-60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, color = 'black', linewidth = 0) +
  annotate('rect',
           xmin = (index.3[1]*10)-60, xmax = (rev(index.3)[1]*10)-60,
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, color = 'black', linewidth = 0) +
  annotate('text', label = 'I', x = -35, y = 0.2, color = 'black', size = font.size-12) +
  annotate('text', label = 'II', x = 50, y = 0.2, color = 'black', size = font.size-12) +
  annotate('text', label = 'III', x = 200, y = 0.2, color = 'black', size = font.size-12) +
  stat_summary(data = df.hpca.sites,
               aes(x = time - 60, y = int, color=site_type, group = site_type),
               fun = median,
               geom = 'line', size = 0.5) +
  stat_summary(data = df.hpca.sites,
               aes(x = time - 60, y = int, color = site_type, group = site_type),
               fun = median,
               geom = 'point', size = 1) +
  stat_summary(data = df.hpca.sites,
               aes(x = time - 60, y = int,
                   color = site_type, fill = site_type, group = site_type),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'ribbon', size = 0, alpha = .25) +
  scale_color_manual(values = c('Spine' = 'red', 'Shaft' = 'green4')) +
  scale_fill_manual(values = c('Spine' = 'red', 'Shaft' = 'green4')) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam, face="bold")) +
  scale_x_continuous(breaks = seq(-60, 230, 20)) +
  labs(caption = 'n = 3/9/81 (cultures/cells/ROIs)',
       x = 'Time, s',
       y = expression(ΔF/F[0]))


# # slope
# df.hpca.sites.med <- df.hpca.sites %>%
#   group_by(site_type, index) %>%
#   mutate(int_med = median(int)) %>%
#   select(-roi_id, -cell_id, -int) %>%
#   distinct() %>%
#   ungroup()
# 
# o1.h = beast(df.hpca.sites.med$int_med[df.hpca.sites.med$site_type == 'Spine'],
#            season = 'none', method='bic')
# df.hpca.sites.beast.spine <- data.frame('Time' = seq(0, length(o1.h$trend$slpSgnPosPr)-1) * 10,
#                                         'Pos.' = o1.h$trend$slpSgnPosPr,
#                                         'Zero' = o1.h$trend$slpSgnZeroPr,
#                                         'Neg.' = rep(1, length(o1.h$trend$slpSgnPosPr)) - (o1.h$trend$slpSgnPosPr+o1.h$trend$slpSgnZeroPr)) %>%
#   pivot_longer(cols = 'Pos.':'Neg.', names_to = 'SlopeSign', values_to = 'p') %>%
#   mutate(SlopeSign = factor(SlopeSign, c('Pos.', 'Zero', 'Neg.'), ordered = TRUE))
# ggplot() +
#   geom_area(data = df.hpca.sites.beast.spine,
#             aes(x = Time, y = p, fill = SlopeSign), alpha = .6) +
#   annotate('text', label = 'NMDA app.', x = 90, y = 1.03, color = 'red') +
#   annotate('rect', xmin = 60, xmax = 120, ymin = 0, ymax = 1.06,
#            alpha = 0.2, fill = 'red') +
#   geom_hline(yintercept = 0.59, lty = 2) +
#   scale_fill_manual(values = c('Pos.' = '#ff4747',
#                                'Zero' = '#49cf36',
#                                'Neg.' =  '#3636cf')) +
#   scale_x_continuous(limits = c(0,280), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(-0.002,1.06), expand = c(0, 0)) +
#   theme_classic() +
#   labs(title = 'Probability of slope sign for HPCA insertions, spines sites',
#        fill = 'Slope sing',
#        x = 'Time, s',
#        y = 'Probability')  
# 
# o2.h = beast(df.hpca.sites.med$int_med[df.hpca.sites.med$site_type == 'Shaft'],
#              season = 'none', method='bic')
# df.hpca.sites.beast.shaft <- data.frame('Time' = seq(0, length(o2.h$trend$slpSgnPosPr)-1) * 10,
#                                         'Pos.' = o2.h$trend$slpSgnPosPr,
#                                         'Zero' = o2.h$trend$slpSgnZeroPr,
#                                         'Neg.' = rep(1, length(o2.h$trend$slpSgnPosPr)) - (o2.h$trend$slpSgnPosPr+o2.h$trend$slpSgnZeroPr)) %>%
#   pivot_longer(cols = 'Pos.':'Neg.', names_to = 'SlopeSign', values_to = 'p') %>%
#   mutate(SlopeSign = factor(SlopeSign, c('Pos.', 'Zero', 'Neg.'), ordered = TRUE))
# ggplot() +
#   geom_area(data = df.hpca.sites.beast.shaft,
#             aes(x = Time, y = p, fill = SlopeSign), alpha = .6) +
#   annotate('text', label = 'NMDA app.', x = 90, y = 1.03, color = 'red') +
#   annotate('rect', xmin = 60, xmax = 120, ymin = 0, ymax = 1.06,
#            alpha = 0.2, fill = 'red') +
#   geom_hline(yintercept = 0.0, lty = 2) +
#   scale_fill_manual(values = c('Pos.' = '#ff4747',
#                                'Zero' = '#49cf36',
#                                'Neg.' =  '#3636cf')) +
#   scale_x_continuous(limits = c(0,280), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(-0.002,1.06), expand = c(0, 0)) +
#   theme_classic() +
#   labs(title = 'Probability of slope sign for HPCA insertions, shaft sites',
#        fill = 'Slope sing',
#        x = 'Time, s',
#        y = 'Probability')  


# time intervals
df.hpca.sites.avg <- df.hpca.sites %>%
  select(-time) %>%
  group_by(roi_id, site_type) %>%
  mutate(time_interval = case_when(index %in% base.indexes ~ 'I',
                                   index %in% mid.indexes.h ~ 'II',
                                   index %in% end.indexes ~ 'III',
                                   .default = 'out')) %>%
  filter(time_interval != 'out') %>%
  droplevels() %>%
  ungroup() %>%
  mutate(time_interval = factor(time_interval, c('I', 'II', 'III'), ordered = TRUE)) %>%
  group_by(roi_id, site_type, time_interval) %>%
  mutate(int_interval = median(int)) %>%
  select(-index, -int) %>%
  ungroup() %>%
  distinct()

# base vs end
df.base.hpca.zero.stat <- df.hpca.sites.avg %>%
  filter(time_interval == 'III') %>%
  select(-time_interval) %>%
  group_by(site_type) %>%
  wilcox_test(int_interval ~ 1, mu = 0) %>%
  add_significance() %>%
  mutate(y.position = c(0.01,0.01), group2 = c(1,1))

plot.box.zero.hpca <- ggplot() +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot(data = df.hpca.sites.avg %>% filter(time_interval == 'I'),
               aes(x = time_interval,
                   y = int_interval,
                   fill = site_type),
               alpha = box.alpha) +
  geom_point(data = df.hpca.sites.avg %>% filter(time_interval == 'I'),
             aes(x = time_interval,
                 y = int_interval),
             size=2, shape=21) +
  stat_pvalue_manual(df.base.hpca.zero.stat, label = 'p.signif',
                     hide.ns = TRUE, remove.bracket = TRUE, size = font.size - 12) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_fill_manual(values = c('Spine' = '#f54040', 'Shaft' = '#4da50b')) +
  scale_y_continuous(limits = c(-0.045, 0.015), breaks = seq(-2, 2, 0.02)) + 
  facet_wrap(~site_type) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam, face = 'bold'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')) +
  labs(x = 'Time interval',
       y = expression(ΔF/F[0]))

# time vs int
df.time.interval.hpca.stat <- df.hpca.sites.avg %>%
  group_by(site_type) %>%
  pairwise_wilcox_test(int_interval ~ time_interval, p.adjust.method = 'BH') %>%
  add_significance() %>%
  add_xy_position(fun = "max") 

plot.box.ti.hpca <- ggplot() +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot(data = df.hpca.sites.avg,
               aes(x = time_interval,
                   y = int_interval,
                   fill = site_type), 
               alpha = box.alpha) +
  geom_point(data = df.hpca.sites.avg,
             aes(x = time_interval,
                 y = int_interval),
             size=2, shape=21) +
  geom_line(data = df.hpca.sites.avg,
            aes(x = time_interval,
                y = int_interval,
                group = roi_id),
            size = .25, lty = 2, alpha = .5) +
  stat_pvalue_manual(df.time.interval.hpca.stat, label = 'p.adj.signif',
                     hide.ns = TRUE, size = font.size - 12) +
  scale_fill_manual(values = c('Spine' = '#f54040', 'Shaft' = '#4da50b')) +
  facet_wrap(~site_type) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam, face="bold")) +
  labs(x = 'Time interval',
       y = expression(ΔF/F[0]))

# site vs int
df.site.type.hpca.stat <- df.hpca.sites.avg %>%
  group_by(time_interval) %>%
  wilcox_test(int_interval ~ site_type) %>%
  add_significance() %>%
  add_xy_position(fun = "max") 

plot.box.si.hpca <- ggplot() +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot(data = df.hpca.sites.avg,
               aes(x = site_type,
                   y = int_interval,
                   fill = site_type),
               alpha = box.alpha) +
  geom_point(data = df.hpca.sites.avg,
             aes(x = site_type,
                 y = int_interval),
             size=2, shape=21) +
  stat_pvalue_manual(df.site.type.hpca.stat, label = 'p.signif',
                     hide.ns = TRUE, size = font.size - 12) +
  scale_fill_manual(values = c('Spine' = '#f54040', 'Shaft' = '#4da50b')) +
  facet_wrap(~time_interval) +
  theme_classic() +
  theme(legend.position = "none",
        text=element_text(size = font.size, family = font.fam, face="bold"),
        axis.text.x = element_text(angle = -90, hjust = 1)) +
  labs(x = 'Site type',
       y = expression(ΔF/F[0]))

##### DF vs Eapp #####

# shaft_rois <- unique(df.fret.sites$roi_id[df.fret.sites$site_type == 'Shaft'])
df.h.vs.f <- df.mask %>%
  filter() %>%
  filter(mask == 'fret',
         index <= 12 & index >= 5,
         (channel == 'ch0' & int_val == 'df') | (channel == 'Eapp' & int_val == 'abs')) %>%
  mutate(roi = as.factor(roi),
         cell_id = id) %>%
  droplevels() %>%
  unite('roi_id', id:roi, sep = '-') %>%
  mutate(roi_id = as.factor(roi_id),
         site_type = as.factor(if_else(roi_id %in% spine_rois, 'Spine', 'Shaft')),
         rel_time = time - 60) %>%
  select(channel, roi_id, int, index, cell_id, site_type, rel_time) %>%
  group_by(roi_id, index) %>%
  pivot_wider(names_from = channel, values_from = int) %>%
  ungroup()


plot.hf <- ggplot(data = df.h.vs.f,
       aes(x = ch0, y = Eapp)) +
  # stat_ellipse(geom="polygon", aes(fill = site_type),
  #              type = "norm", alpha = .3, level = 0.95) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_mark_hull(aes(fill = site_type, label = site_type),
                    size = 0, expand = 0.0075) + 
  # geom_line(aes(color = rel_time, group = roi_id),
  #           alpha = .5) +
  geom_point(aes(color = rel_time, shape = site_type), alpha = .6, size = 3) +
  scale_fill_manual(values = c('Spine' = '#f54040', 'Shaft' = '#4da50b')) + 
  scale_color_gradient2(low = "green2", mid = 'yellow', high = "red4") +
  theme_classic() +
  theme(text=element_text(size = font.size - 3, family = font.fam, face="bold")) +
  labs(color = 'Time',
       x = expression(ΔF/F[0]),
       y = expression(E[app])) +
  guides(fill = 'none',
         shape = 'none')
plot.hf

draw.hf <- ggdraw(plot.hf) +
  draw_plot_label(c("F"),
                  c(0),
                  c(1),
                  size = font.size + 3)
draw.hf
save_plot('0_plot_h_vs_f.png', draw.hf, base_width = 12, base_height = 4)



##### PLOTS #####
# fret hist
f.hist.b <- plot_hist(df.h = df.fret.base,
                      scale = 0.5,
                      df.g = df.base_mixmdl,
                      df.b = df.gmm.optim  %>% filter(group == 'fret_base'),
                      mm = base_mixmgl,
                      bic.title = 'E')

f.hist.b
save_plot('0_plot_fret_hist.png', f.hist.b, base_width = 7, base_height = 5)

# fret
draw.fret.abs <- ggdraw(plot.fret.abs) +
  draw_plot(plot.box.zero.abs, x = 0.1, y = 0.58, width = 0.25, height = 0.3) +
  draw_plot_label(c("Ha"),
                  c(0),
                  c(0.91),
                  size = font.size + 3)
draw.fret.df <- ggdraw(plot.fret.df) +
  draw_plot_label(c("Ia"),
                  c(0),
                  c(0.91),
                  size = font.size + 3)
# draw.box.zero <- ggdraw(plot.box.zero.abs) +
#   draw_plot_label(c("E"),
#                   c(0),
#                   c(1),
#                   size = font.size + 3)
draw.box.si.abs <- ggdraw(plot.box.si.abs) +
  draw_plot_label(c("Hc"),
                  c(0),
                  c(1),
                  size = font.size + 3)
draw.box.ti.abs <- ggdraw(plot.box.ti.abs) +
  draw_plot_label(c("Hb"),
                  c(0),
                  c(1),
                  size = font.size + 3)
draw.box.si.df <- ggdraw(plot.box.si.df) +
  draw_plot_label(c("Ic"),
                  c(0),
                  c(1),
                  size = font.size + 3)
draw.box.ti.df <- ggdraw(plot.box.ti.df) +
  draw_plot_label(c("Ib"),
                  c(0),
                  c(1),
                  size = font.size + 3)

fret.abs <- plot_grid(draw.fret.abs, draw.box.ti.abs, draw.box.si.abs,
          rel_widths = c(1,0.5,0.5), nrow=1)
save_plot('0_plot_fret_abs.png', fret.abs, base_width = 18, base_height = 4)

fret.df <- plot_grid(draw.fret.df, draw.box.ti.df, draw.box.si.df,
          rel_widths = c(1,0.5,0.5), nrow=1)
save_plot('0_plot_fret_df.png', fret.df, base_width = 18, base_height = 4)

# hpca
draw.hpca.df <- ggdraw(plot.hpca.df) +
  draw_plot(plot.box.zero.hpca, x = 0.65, y = 0.6, width = 0.25, height = 0.3) +
  draw_plot_label(c("Ga"),
                  c(0),
                  c(0.91),
                  size = font.size + 3)
draw.box.si.hpca <- ggdraw(plot.box.si.hpca) +
  draw_plot_label(c("Gb"),
                  c(0),
                  c(1),
                  size = font.size + 3)
draw.box.ti.hpca <- ggdraw(plot.box.ti.hpca) +
  draw_plot_label(c("Gc"),
                  c(0),
                  c(1),
                  size = font.size + 3)

hpca.df <- plot_grid(draw.hpca.df, draw.box.ti.hpca, draw.box.si.hpca,
          rel_widths = c(1,0.5,0.5), nrow=1)
save_plot('0_plot_hpca_df.png', hpca.df, base_width = 18, base_height = 4)
