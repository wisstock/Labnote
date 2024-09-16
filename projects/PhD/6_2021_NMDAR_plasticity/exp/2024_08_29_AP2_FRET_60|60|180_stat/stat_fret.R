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
  return(list(data.frame(val = input_dens$x,
                    raw = input_dens$y,
                    comp1 = comp1,
                    comp2 = comp2,
                    comp_comb = comp1+comp2),
              mixmdl))
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
base.indexes <- c(0,1,2,3,4,5,6)
mid.indexes <- c(12,13,14,15)
end.indexes <- c(26,27,28)

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
           alpha = 0.075, color = 'black') +
  annotate('rect',
           xmin = mid.indexes[1], xmax = rev(mid.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black') +
  annotate('rect',
           xmin = end.indexes[1], xmax = rev(end.indexes)[1],
           ymin = -Inf, ymax = Inf,
           alpha = 0.075, color = 'black') +
  labs(title = 'FRET profiles in individual ROIs',
       caption = 'Red rect - NMDA app., black rect - pre and post app. frames',
       y = expression(E[app]),
       x = 'Frame idx') +
  scale_fill_manual(values = rainbow(length(df.fret.plot$roi_id))) +
  theme_classic() +
  theme(legend.position="none") +
  facet_wrap(~cell_id)

# ggsave('fret_prof.png', p1, dpi = 300)


# base
df.fret.base <- df.mask %>%
  filter(channel == 'Eapp',
         index %in% base.indexes) %>%
  select(id, roi, int, mask) %>%
  droplevels() %>%
  group_by(id, roi, mask) %>%
  mutate(int = mean(int),
         roi = as.factor(roi),
         cell_id = id) %>%
  ungroup() %>%
  distinct() %>%
  unite('roi_id', id:roi, sep = '-')

bm <- calc.mixmdl(df.fret.base$int[df.fret.base$mask == selected.mask]) 
df.base_mixmdl <- bm[[1]]
base_mixmgl <- bm[[2]]
remove(bm)
# ggplot(data = df.pre_mixmdl, aes(x = val)) +
#   geom_line(aes(y = comp1), color = 'red', alpha = .5) +
#   geom_line(aes(y = comp2), color = 'blue', alpha = .5) +
#   geom_line(aes(y = comp_comb), lty = 2) +
#   geom_line(aes(y = raw))

plot.base <- ggplot() +
  geom_dotplot(data = df.fret.base %>% filter(mask == selected.mask),
               aes(x = int, fill = roi_id),
               stackgroups = TRUE, binwidth = b.width,
               binpositions = "all", method = "histodot") +
  geom_line(data = df.base_mixmdl,
            aes(x = val, y = comp1 / comp.scale.factor),
            color = 'red') +
  geom_line(data = df.base_mixmdl,
            aes(x = val, y = comp2 / comp.scale.factor),
            color = 'blue') +
  geom_line(data = df.base_mixmdl,
            aes(x = val, y = comp_comb / comp.scale.factor),
            lty = 2) +
  geom_vline(xintercept = base_mixmgl$mu, lty = 3) +
  scale_fill_manual(values = rainbow(length(df.fret.base$roi_id[df.fret.base$mask == selected.mask]))) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank()) +
  labs(title = 'Base app. average FRET in individual ROIs',
       x = expression(E[app])) +
  xlim(-0.005,0.11)

# ggsave('fret_pre.png', p2, dpi = 300)

# mid
df.fret.mid <- df.mask %>%
  filter(channel == 'Eapp',
         index %in% mid.indexes) %>%
  select(id, roi, int, mask) %>%
  droplevels() %>%
  group_by(id, roi, mask) %>%
  mutate(int = mean(int),
         roi = as.factor(roi),
         cell_id = id) %>%
  ungroup() %>%
  distinct() %>%
  unite('roi_id', id:roi, sep = '-')

mm <- calc.mixmdl(df.fret.mid$int[df.fret.mid$mask == selected.mask]) 
df.mid_mixmdl <- mm[[1]]
mid_mixmgl <- mm[[2]]
remove(mm)

plot.mid <- ggplot() +
  geom_dotplot(data = df.fret.mid %>% filter(mask == selected.mask),
               aes(x = int, fill = roi_id),
               stackgroups = TRUE, binwidth = b.width,
               binpositions = "all", method = "histodot") +
  geom_line(data = df.mid_mixmdl,
            aes(x = val, y = comp1 / comp.scale.factor),
            color = 'red') +
  geom_line(data = df.mid_mixmdl,
            aes(x = val, y = comp2 / comp.scale.factor),
            color = 'blue') +
  geom_line(data = df.mid_mixmdl,
            aes(x = val, y = comp_comb / comp.scale.factor),
            lty = 2) +
  geom_vline(xintercept = mid_mixmgl$mu, lty = 3) +
  scale_fill_manual(values = rainbow(length(df.fret.mid$roi_id[df.fret.mid$mask == selected.mask]))) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank()) +
  labs(title = 'Mid app. average FRET in individual ROIs',
       x = expression(E[app])) +
  xlim(-0.005,0.11)

# ggsave('fret_post.png', p3, dpi = 300)


# end
df.fret.end <- df.mask %>%
  filter(channel == 'Eapp',
         index %in% end.indexes) %>%
  select(id, roi, int, mask) %>%
  droplevels() %>%
  group_by(id, roi, mask) %>%
  mutate(int = mean(int),
         roi = as.factor(roi),
         cell_id = id) %>%
  ungroup() %>%
  distinct() %>%
  unite('roi_id', id:roi, sep = '-')

em <- calc.mixmdl(df.fret.end$int[df.fret.end$mask == selected.mask]) 
df.end_mixmdl <- em[[1]]
end_mixmgl <- em[[2]]
remove(em)

plot.end <- ggplot() +
  geom_dotplot(data = df.fret.end %>% filter(mask == selected.mask),
               aes(x = int, fill = roi_id),
               stackgroups = TRUE, binwidth = b.width,
               binpositions = "all", method = "histodot") +
  geom_line(data = df.end_mixmdl,
            aes(x = val, y = comp1 / comp.scale.factor),
            color = 'red') +
  geom_line(data = df.end_mixmdl,
            aes(x = val, y = comp2 / comp.scale.factor),
            color = 'blue') +
  geom_line(data = df.end_mixmdl,
            aes(x = val, y = comp_comb / comp.scale.factor),
            lty = 2) +
  geom_vline(xintercept = end_mixmgl$mu, lty = 3) +
  scale_fill_manual(values = rainbow(length(df.fret.end$roi_id[df.fret.end$mask == selected.mask]))) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank()) +
  labs(title = 'End app. average FRET in individual ROIs',
       x = expression(E[app])) +
  xlim(-0.005,0.11)


# combined plot
plot_grid(plot.base, plot.mid, plot.end,
          ncol = 1)  # labels = "AUTO",


##### TIMEPOINTS BOX #####
sites.threshold <- 0.03

df.fret.sites <- df.mask %>%
  filter(channel == 'Eapp',
         int_val == 'abs',
         mask == 'fret') %>%
  mutate(roi = as.factor(roi),
         cell_id = id) %>%
  droplevels() %>%
  unite('roi_id', id:roi, sep = '-') %>%
  group_by(roi_id) %>%
  mutate(wash_tail_int = mean(int[index %in% seq(16,28)]),
         site_type = as.factor(ifelse(wash_tail_int > sites.threshold, 'spine', 'shaft'))) %>%
  ungroup() %>%
  select(-time, -dist, -dist_group, -int_val, -channel, -mask)




# df.to.box <- df.fret.sites  %>%  # by time points
#   filter(index %in% c(3, 15, 24)) %>%
#   droplevels() %>%
#   select(-roi_id, -wash_tail_int) %>%
#   mutate(index = as.factor(index)) %>%
#   group_by(cell_id, index, site_type) %>%
#   mutate(int_med = median(int), int_iqr = IQR(int)) %>%
#   select(-int) %>%
#   ungroup() %>%
#   distinct()

df.to.box <- df.fret.sites  %>%  # by time points
 filter(index %in% c(3, 15, 24)) %>%
 droplevels() %>%
 select(-roi_id, -wash_tail_int) %>%
 mutate(index = as.factor(index)) %>%
 group_by(cell_id, index, site_type) %>%
 mutate(int_med = median(int), int_iqr = IQR(int)) %>%
 select(-int) %>%
 ungroup() %>%
 distinct()


# IS ZERO
wilcox.test(df.to.box$int_med[df.to.box$index == 3 & df.to.box$site_type == 'spine'], mu = 0)
wilcox.test(df.to.box$int_med[df.to.box$index == 3 & df.to.box$site_type == 'shaft'], mu = 0)

df.base.zero.stat <- df.to.box %>%
  filter(index == 3) %>%
  group_by(site_type) %>%
  wilcox_test(int_med~1, mu = 0) %>%
  add_significance() %>%
  mutate(y.position = c(0.045,0.02), group2 = c(0,0))

ggplot(df.to.box %>% filter(index == 3),
       aes(y = int_med)) +
  geom_boxplot() +
  stat_pvalue_manual(df.base.zero.stat, label = 'p.signif',
                     hide.ns = TRUE, remove.bracket = TRUE, label.size = 7) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = expression(E[app])) +
  facet_wrap(~site_type)
  


# I vs T
df.to.box.stat.it <- df.to.box %>%
  select(-cell_id, -int_iqr) %>%
  group_by(site_type) %>%
  pairwise_wilcox_test(int_med ~ index, p.adjust.method = "BH") %>%
  add_xy_position(fun = "mean_sd")

ggplot(data = df.to.box,  #  %>% filter(mask == 'fret')
       aes(x = index, y = int_med)) +
  geom_boxplot(aes(fill = index)) +
  stat_pvalue_manual(df.to.box.stat.it, label = 'p.adj.signif',
                     hide.ns = TRUE) +
  facet_wrap(~site_type, ncol = 2) 

# I vs M
df.to.box.stat.im <- df.to.box %>%
  select(-cell_id, -int_iqr) %>%
  group_by(index) %>%
  wilcox_test(int_med ~ site_type) %>%
  add_significance() %>%
  add_xy_position(fun = "mean_sd")  

ggplot(data = df.to.box,
       aes(x = site_type, y = int_med)) +
  geom_boxplot(aes(fill = site_type)) +
  stat_pvalue_manual(df.to.box.stat.im, label = 'p.signif',
                     hide.ns = TRUE) +
  facet_wrap(~index, ncol = 4)