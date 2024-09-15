# NMDA ionophoresis, FRET between AP2B1-EYFP and HPCA(WT)-ECFP
# Copyright © 2024 Borys Olifirov

require(stringr)
require(dplyr)
require(tidyr)
require(purrr)
require(mixtools)
require(rstatix)
require(ggplot2)
require(ggpubr)
require(cowplot)
require(ggsci)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2024_08_29_AP2_FRET_60|60|180_stat')

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
  filter(index <= 29)  # , , id != '24_05_16_09' , id != '24_06_5_04'

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


##### LOW vs HIGH FRET #####
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
  return(data.frame(val = input_dens$x,
                    raw = input_dens$y,
                    comp1 = comp1,
                    comp2 = comp2,
                    comp_comb = comp1+comp2))
}

calc.mixmdl3 <- function(input_vector) {
  mixmdl <- normalmixEM(input_vector, k = 3)
  input_dens <- density(input_vector)
  comp1 <- dnorm(x = input_dens$x,
                 mean = mixmdl$mu[1],
                 sd = mixmdl$sigma[1]) * mixmdl$lambda[1]
  comp2 <- dnorm(x = input_dens$x,
                 mean = mixmdl$mu[2],
                 sd = mixmdl$sigma[2]) * mixmdl$lambda[2]
  comp3 <- dnorm(x = input_dens$x,
                 mean = mixmdl$mu[3],
                 sd = mixmdl$sigma[3]) * mixmdl$lambda[3]
  return(data.frame(val = input_dens$x,
                    raw = input_dens$y,
                    comp1 = comp1,
                    comp2 = comp2,
                    comp3 = comp3,
                    comp_comb = comp1+comp2+comp3))
}

selected.mask <- 'fret'
comp.scale.factor <- 50
b.width <- 0.003
pre.indexes <- c(0,1,2,3,4,5,6)
post.indexes <- c(12,13,14,15)

# selected frames
df.fret.plot <- df.mask %>%
  filter(channel == 'Eapp',
         int_val == 'abs',
         mask == selected.mask) %>%
  mutate(roi = as.factor(roi),
         cell_id = id) %>%
  droplevels() %>%
  unite('roi_id', id:roi, sep = '-')

p1 <- ggplot(data = df.fret.plot,
       aes(x = index, y = int, color = roi_id, group = roi_id)) +
  stat_summary(fun = median,
               geom = 'line', size = .5) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(data = df.fret.plot,
               aes(x = index, y = int, group = channel),
               color = 'black',
               fun = median,
               geom = 'point', size = 0.75) +
  stat_summary(data = df.fret.plot,
               aes(x = index, y = int, group = channel),
               color = 'black',
               fun = median,
               geom = 'line', size = 0.3) +
  stat_summary(data = df.fret.plot,
               aes(x = index, y = int, group = channel),
               color = 'black',
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'errorbar', size = 0.15, width = 0.75) +
  annotate('rect', xmin = 6, xmax = 12, ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = 'red') +
  annotate('rect',
           xmin = pre.indexes[1], xmax = rev(pre.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black') +
  annotate('rect',
           xmin = post.indexes[1], xmax = rev(post.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black') +
  labs(title = 'FRET profiles in individual ROIs',
       caption = 'Red rect - NMDA app., black rect - pre and post app. frames',
       y = expression(E[app]),
       x = 'Frame idx') +
  scale_fill_manual(values = rainbow(length(df.fret.plot$roi_id))) +
  theme_classic() +
  theme(legend.position="none") +
  facet_wrap(~id)

ggsave('fret_prof.png', p1, dpi = 300)


# pre
df.fret.pre <- df.mask %>%
  filter(channel == 'Eapp',
         index %in% pre.indexes) %>%
  select(id, roi, int, mask) %>%
  droplevels() %>%
  group_by(id, roi, mask) %>%
  mutate(int = mean(int),
         roi = as.factor(roi),
         cell_id = id) %>%
  ungroup() %>%
  distinct() %>%
  unite('roi_id', id:roi, sep = '-')

df.pre_mixmdl <- calc.mixmdl(df.fret.pre$int[df.fret.pre$mask == selected.mask]) 
# ggplot(data = df.pre_mixmdl, aes(x = val)) +
#   geom_line(aes(y = comp1), color = 'red', alpha = .5) +
#   geom_line(aes(y = comp2), color = 'blue', alpha = .5) +
#   geom_line(aes(y = comp_comb), lty = 2) +
#   geom_line(aes(y = raw))

p2 <- ggplot() +
  geom_dotplot(data = df.fret.pre %>% filter(mask == selected.mask),
               aes(x = int, fill = roi_id),
               stackgroups = TRUE, binwidth = b.width,
               binpositions = "all", method = "histodot") +
  geom_line(data = df.pre_mixmdl,
            aes(x = val, y = comp1 / comp.scale.factor),
            color = 'red') +
  geom_line(data = df.pre_mixmdl,
            aes(x = val, y = comp2 / comp.scale.factor),
            color = 'blue') +
  geom_line(data = df.pre_mixmdl,
            aes(x = val, y = comp_comb / comp.scale.factor),
            lty = 2) +
  scale_fill_manual(values = rainbow(length(df.fret.post$roi_id[df.fret.post$mask == selected.mask]))) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank()) +
  labs(title = 'Pre app. average FRET in individual ROIs',
       x = expression(E[app])) +
  xlim(-0.005,0.13)

ggsave('fret_pre.png', p2, dpi = 300)

# post
df.fret.post <- df.mask %>%
  filter(channel == 'Eapp',
         index %in% post.indexes) %>%
  select(id, roi, int, mask) %>%
  droplevels() %>%
  group_by(id, roi, mask) %>%
  mutate(int = mean(int),
         roi = as.factor(roi),
         cell_id = id) %>%
  ungroup() %>%
  distinct() %>%
  unite('roi_id', id:roi, sep = '-')

df.post_mixmdl <- calc.mixmdl(df.fret.post$int[df.fret.post$mask == selected.mask]) 

p3 <- ggplot() +
  geom_dotplot(data = df.fret.post %>% filter(mask == selected.mask),
               aes(x = int, fill = roi_id),
               stackgroups = TRUE, binwidth = b.width,
               binpositions = "all", method = "histodot") +
  geom_line(data = df.post_mixmdl,
            aes(x = val, y = comp1 / comp.scale.factor),
            color = 'red') +
  geom_line(data = df.post_mixmdl,
            aes(x = val, y = comp2 / comp.scale.factor),
            color = 'blue') +
  geom_line(data = df.post_mixmdl,
            aes(x = val, y = comp_comb / comp.scale.factor),
            lty = 2) +
  scale_fill_manual(values = rainbow(length(df.fret.post$roi_id[df.fret.post$mask == selected.mask]))) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank()) +
  labs(title = 'Post app. average FRET in individual ROIs',
       x = expression(E[app])) +
  xlim(-0.005,0.13)

ggsave('fret_post.png', p3, dpi = 300)


##### TIMEPOINTS BOX #####
df.to.box <- df.mask %>%
  filter(index %in% c(2, 9, 13, 25),
         channel == 'Eapp',
         dist_group != 'min',
         id != '24_05_16_04') %>%
  droplevels() %>%
  group_by(id, channel, index, mask) %>%
  mutate(int_med = median(int), int_iqr = IQR(int), idx = as.factor(index)) %>%
  ungroup() %>%
  select(-time, -roi, -int, -dist, -dist_group, -channel, -index) %>%
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

