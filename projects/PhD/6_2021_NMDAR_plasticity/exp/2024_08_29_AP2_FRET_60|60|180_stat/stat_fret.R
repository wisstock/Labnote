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


calc.mixmdl.optim <- function(input_vector) {
  comp.vals <- seq(2,10)
  loglik.vector <- c()
  for (comp in comp.vals) {
    mixmdl <- normalmixEM(input_vector, k = comp)
    loglik.vector <- append(loglik.vector, mixmdl$loglik)
  }
  return(loglik.vector)
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
base.indexes <- seq(0,5)
mid.indexes.f <- seq(10,12)
mid.indexes.h <- seq(7,9)
end.indexes <- seq(24,28)
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
  labs(title = 'FRET profiles in individual ROIs',
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
  labs(title = 'HPCA insertions profiles in individual ROIs',
       caption = 'Red rect - NMDA app., black rect - time intervals for analysis',
       y = expression(ΔF/F[0]),
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

df.hpca.base <- df.mask %>%
  filter(channel == 'ch0', int_val == 'df',
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


comp.scale.factor <- 3
b.width <- 0.003

plot(calc.mixmdl.optim(df.fret.base$int[df.fret.base$mask == selected.mask]), type = 'line')
plot(calc.mixmdl.optim(df.hpca.base$int[df.hpca.base$mask == selected.mask]), type = 'line')

bm <- calc.mixmdl(df.fret.base$int[df.fret.base$mask == selected.mask]) 
df.base_mixmdl <- bm[[1]]
base_mixmgl <- bm[[2]]
remove(bm)

bm.h <- calc.mixmdl(df.hpca.base$int[df.hpca.base$mask == selected.mask]) 
df.base_mixmdl.h <- bm.h[[1]]
base_mixmgl_h <- bm.h[[2]]
remove(bm.h)
# ggplot(data = df.pre_mixmdl, aes(x = val)) +
#   geom_line(aes(y = comp1), color = 'red', alpha = .5) +
#   geom_line(aes(y = comp2), color = 'blue', alpha = .5) +
#   geom_line(aes(y = comp_comb), lty = 2) +
#   geom_line(aes(y = raw))

plot.base <- ggplot() +
  geom_histogram(data = df.fret.base %>% filter(mask == selected.mask),
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


plot.base.h <- ggplot() +
  geom_histogram(data = df.hpca.base %>% filter(mask == selected.mask),
                 aes(x = int, fill = roi_id),
                 stackgroups = TRUE, binwidth = b.width,
                 binpositions = "all", method = "histodot") +
  geom_line(data = df.base_mixmdl.h,
            aes(x = val, y = comp1 / comp.scale.factor),
            color = 'red') +
  geom_line(data = df.base_mixmdl.h,
            aes(x = val, y = comp2 / comp.scale.factor),
            color = 'blue') +
  geom_line(data = df.base_mixmdl.h,
            aes(x = val, y = comp_comb / comp.scale.factor),
            lty = 2) +
  geom_vline(xintercept = base_mixmgl_h$mu, lty = 3) +
  scale_fill_manual(values = rainbow(length(df.hpca.base$roi_id[df.hpca.base$mask == selected.mask]))) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank()) +
  labs(title = 'Base app. average HPCA insertions in individual ROIs',
       x = expression(ΔF/F[0])) +
  xlim(-0.5,1)


# ggsave('fret_pre.png', p2, dpi = 300)

# mid
df.fret.mid <- df.mask %>%
  filter(channel == 'Eapp',
         index %in% mid.indexes.f) %>%
  select(id, roi, int, mask) %>%
  droplevels() %>%
  group_by(id, roi, mask) %>%
  mutate(int = mean(int),
         roi = as.factor(roi),
         cell_id = id) %>%
  ungroup() %>%
  distinct() %>%
  unite('roi_id', id:roi, sep = '-')

df.hpca.mid <- df.mask %>%
  filter(channel == 'ch0', int_val == 'df',
         index %in% mid.indexes.h) %>%
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

mm.h <- calc.mixmdl(df.hpca.mid$int[df.hpca.mid$mask == selected.mask]) 
df.mid_mixmdl.h <- mm.h[[1]]
mid_mixmgl_h <- mm.h[[2]]
remove(mm.h)

plot.mid <- ggplot() +
  geom_histogram(data = df.fret.mid %>% filter(mask == selected.mask),
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

plot.mid.h <- ggplot() +
  geom_histogram(data = df.hpca.mid %>% filter(mask == selected.mask),
                 aes(x = int, fill = roi_id),
                 stackgroups = TRUE, binwidth = 0.03,
                 binpositions = "all", method = "histodot") +
  geom_line(data = df.mid_mixmdl.h,
            aes(x = val, y = comp1 / 0.3),
            color = 'red') +
  geom_line(data = df.mid_mixmdl.h,
            aes(x = val, y = comp2 / 0.3),
            color = 'blue') +
  geom_line(data = df.mid_mixmdl.h,
            aes(x = val, y = comp_comb / 0.3),
            lty = 2) +
  geom_vline(xintercept = mid_mixmgl_h$mu, lty = 3) +
  scale_fill_manual(values = rainbow(length(df.hpca.mid$roi_id[df.hpca.mid$mask == selected.mask]))) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank()) +
  labs(title = 'Mid app. average HPCA insertions in individual ROIs',
       x = expression(ΔF/F[0])) +
  xlim(-0.5,1)

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

df.hpca.end <- df.mask %>%
  filter(channel == 'ch0', int_val == 'df',
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

em.h <- calc.mixmdl(df.hpca.end$int[df.hpca.end$mask == selected.mask]) 
df.end_mixmdl.h <- em.h[[1]]
end_mixmgl_h <- em.h[[2]]
remove(em.h)

plot.end <- ggplot() +
  geom_histogram(data = df.fret.end %>% filter(mask == selected.mask),
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


plot.end.h <- ggplot() +
  geom_histogram(data = df.hpca.end %>% filter(mask == selected.mask),
                 aes(x = int, fill = roi_id),
                 stackgroups = TRUE, binwidth = 0.03,
                 binpositions = "all", method = "histodot") +
  geom_line(data = df.end_mixmdl.h,
            aes(x = val, y = comp1 / 0.3),
            color = 'red') +
  geom_line(data = df.end_mixmdl.h,
            aes(x = val, y = comp2 / 0.3),
            color = 'blue') +
  geom_line(data = df.end_mixmdl.h,
            aes(x = val, y = comp_comb / 0.3),
            lty = 2) +
  geom_vline(xintercept = end_mixmgl_h$mu, lty = 3) +
  scale_fill_manual(values = rainbow(length(df.hpca.end$roi_id[df.hpca.end$mask == selected.mask]))) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y = element_blank()) +
  labs(title = 'End app. average HPCA insertions in individual ROIs',
       x = expression(ΔF/F[0])) +
  xlim(-0.5,1)


# combined plot
plot_grid(plot.base, plot.base.h,
          plot.mid, plot.mid.h,
          plot.end, plot.end.h,
          ncol = 2)  # labels = "AUTO",


df.for.gmm <- data.frame('fret_base' = df.fret.base$int[df.fret.base$mask == selected.mask],
                         'hpca_base' = df.hpca.base$int[df.hpca.base$mask == selected.mask],
                         'fret_mid' = df.fret.mid$int[df.fret.mid$mask == selected.mask],
                         'hpca_mid' = df.hpca.mid$int[df.hpca.mid$mask == selected.mask],
                         'fret_end' = df.fret.end$int[df.fret.end$mask == selected.mask],
                         'hpca_end' = df.hpca.end$int[df.hpca.end$mask == selected.mask])
write.csv(df.for.gmm, file = 'fret_for_gmm.csv')

plot(density(df.hpca.mid$int[df.hpca.mid$mask == selected.mask]))

##### TIME INTERVALS FRET #####
sites.threshold <- 0.03
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
  mutate(wash_tail_int = mean(int[index %in% seq(16,28)]),
         site_type = as.factor(ifelse(wash_tail_int > sites.threshold, 'Spine', 'Shaft'))) %>%
  ungroup() %>%
  select(-dist, -dist_group, -int_val, -channel, -mask)

# profiles
ggplot() +
  stat_summary(data = df.fret.sites,
               aes(x = time, y = int, color=site_type, group = site_type),
               fun = median,
               geom = 'line', size = 0.5) +
  stat_summary(data = df.fret.sites,
               aes(x = time, y = int, color = site_type, group = site_type),
               fun = median,
               geom = 'point', size = 1) +
  stat_summary(data = df.fret.sites,
               aes(x = time, y = int,
                   color = site_type, fill = site_type, group = site_type),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'ribbon', size = 0, alpha = .25) +
  annotate('text', label = 'NMDA app.', x = 100, y = 0.075, color = 'red') +
  annotate('rect', xmin = 60, xmax = 120, ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = 'red') +
  annotate('text', label = 'Base.', x = 30, y = 0.0675, color = 'black') +
  annotate('rect',
           xmin = index.1[1] * 10, xmax = rev(index.1)[1] * 10,
           ymin = -Inf, ymax = 0.07,
           alpha = 0.075, color = 'black', linewidth = 0) +
  annotate('text', label = 'Mid.', x = 110, y = 0.0675, color = 'black') +
  annotate('rect',
           xmin = index.2[1] * 10, xmax = rev(index.2)[1] * 10,
           ymin = -Inf, ymax = 0.07,
           alpha = 0.075, color = 'black', linewidth = 0) +
  annotate('text', label = 'End', x = 260, y = 0.0675, color = 'black') +
  annotate('rect',
           xmin = index.3[1] * 10, xmax = rev(index.3)[1] * 10,
           ymin = -Inf, ymax = 0.07,
           alpha = 0.075, color = 'black', linewidth = 0) +
  scale_color_manual(values = c('Spine' = '#f54040', 'Shaft' = '#4da50b')) +
  scale_fill_manual(values = c('Spine' = '#f54040', 'Shaft' = '#4da50b')) +
  labs(title = 'FRET in different ROI types',
       color = 'ROI type',
       fill = 'ROI type',
       x = 'Time, s',
       y = expression(E[app])) +
  theme_classic()


# plat det
df.fret.sites.med <- df.fret.sites %>%
  group_by(site_type, index) %>%
  mutate(int_med = median(int)) %>%
  select(-roi_id, -cell_id, -wash_tail_int, -int) %>%
  distinct() %>%
  ungroup()

o1 = beast(df.fret.sites.med$int_med[df.fret.sites.med$site_type == 'Spine'],
           season = 'none', method='bic')
df.fret.sites.beast.spine <- data.frame('Time' = seq(0, length(o1$trend$slpSgnPosPr)-1) * 10,
                                        'Pos.' = o1$trend$slpSgnPosPr,
                                        'Zero' = o1$trend$slpSgnZeroPr,
                                        'Neg.' = rep(1, length(o1$trend$slpSgnPosPr)) - (o1$trend$slpSgnPosPr+o1$trend$slpSgnZeroPr)) %>%
  pivot_longer(cols = 'Pos.':'Neg.', names_to = 'SlopeSign', values_to = 'p') %>%
  mutate(SlopeSign = factor(SlopeSign, c('Pos.', 'Zero', 'Neg.'), ordered = TRUE))
ggplot() +
  geom_area(data = df.fret.sites.beast.spine,
            aes(x = Time, y = p, fill = SlopeSign), alpha = .6) +
  annotate('text', label = 'NMDA app.', x = 90, y = 1.03, color = 'red') +
  annotate('rect', xmin = 60, xmax = 120, ymin = 0, ymax = 1.06,
           alpha = 0.2, fill = 'red') +
  geom_hline(yintercept = 0.82, lty = 2) +
  scale_fill_manual(values = c('Pos.' = '#ff4747',
                               'Zero' = '#49cf36',
                               'Neg.' =  '#3636cf')) +
  scale_x_continuous(limits = c(0,280), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.002,1.06), expand = c(0, 0)) +
  theme_classic() +
  labs(title = 'Probability of slope sign for spines sites',
       fill = 'Slope sing',
       x = 'Time, s',
       y = 'Probability')  

o2 = beast(df.fret.sites.med$int_med[df.fret.sites.med$site_type == 'Shaft'],
           season = 'none', method='bic')
df.fret.sites.beast.shaft <- data.frame('Time' = seq(0, length(o2$trend$slpSgnPosPr)-1) * 10,
                                        'Pos.' = o2$trend$slpSgnPosPr,
                                        'Zero' = o2$trend$slpSgnZeroPr,
                                        'Neg.' = rep(1, length(o2$trend$slpSgnPosPr)) - (o2$trend$slpSgnPosPr+o2$trend$slpSgnZeroPr)) %>%
  pivot_longer(cols = 'Pos.':'Neg.', names_to = 'SlopeSign', values_to = 'p') %>%
  mutate(SlopeSign = factor(SlopeSign, c('Pos.', 'Zero', 'Neg.'), ordered = TRUE))
ggplot() +
  geom_area(data = df.fret.sites.beast.shaft,
            aes(x = Time, y = p, fill = SlopeSign), alpha = .6) +
  annotate('text', label = 'NMDA app.', x = 90, y = 1.03, color = 'red') +
  annotate('rect', xmin = 60, xmax = 120, ymin = 0, ymax = 1.06,
           alpha = 0.2, fill = 'red') +
  geom_hline(yintercept = 1, lty = 2) +
  scale_fill_manual(values = c('Pos.' = '#ff4747',
                               'Zero' = '#49cf36',
                               'Neg.' =  '#3636cf')) +
  scale_x_continuous(limits = c(0,280), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.002,1.06), expand = c(0, 0)) +
  theme_classic() +
  labs(title = 'Probability of slope sign for shaft sites',
       fill = 'Slope sing',
       x = 'Time, s',
       y = 'Probability')  

# time intervals
df.fret.sites.avg <- df.fret.sites %>%
  select(-wash_tail_int, -time) %>%
  group_by(roi_id, site_type) %>%
  mutate(time_interval = case_when(index %in% index.1 ~ 'base',
                                   index %in% index.2 ~ 'mid',
                                   index %in% index.3 ~ 'end',
                                   .default = 'out')) %>%
  filter(time_interval != 'out') %>%
  droplevels() %>%
  ungroup() %>%
  mutate(time_interval = factor(time_interval, c('base', 'mid', 'end'), ordered = TRUE)) %>%
  group_by(roi_id, site_type, time_interval) %>%
  mutate(int_interval = median(int)) %>%
  select(-index, -int) %>%
  ungroup() %>%
  distinct()

# base vs zero
df.base.zero.stat <- df.fret.sites.avg %>%
  filter(time_interval == 'base') %>%
  select(-time_interval) %>%
  group_by(site_type) %>%
  wilcox_test(int_interval ~ 1, mu = 0) %>%
  add_significance() %>%
  mutate(y.position = c(0.082,0.035), group2 = c(1,1))

ggplot() +
  geom_boxplot(data = df.fret.sites.avg %>% filter(time_interval == 'base'),
               aes(x = time_interval,
                   y = int_interval,
                   fill = site_type)) +
  geom_point(data = df.fret.sites.avg %>% filter(time_interval == 'base'),
             aes(x = time_interval,
                 y = int_interval),
             size=2, shape=21) +
  stat_pvalue_manual(df.base.zero.stat, label = 'p.signif',
                     hide.ns = TRUE, remove.bracket = TRUE, label.size = 5) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c('Spine' = '#f54040', 'Shaft' = '#4da50b')) +
  labs(title = 'Initial FRET in different ROI types',
       y = expression(E[app])) +
  facet_wrap(~site_type)

# time vs int
df.time.interval.stat <- df.fret.sites.avg %>%
  group_by(site_type) %>%
  pairwise_wilcox_test(int_interval ~ time_interval, p.adjust.method = 'BH') %>%
  add_significance() %>%
  add_xy_position(fun = "max") 

ggplot() +
  geom_boxplot(data = df.fret.sites.avg,
               aes(x = time_interval,
                   y = int_interval,
                   fill = site_type)) +
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
                     hide.ns = TRUE) +
  theme_classic() +
  theme(legend.position="none") +
  labs(title = 'Time courses of FRET in different ROI types',
       x = 'Time interval',
       y = expression(E[app])) +
  scale_fill_manual(values = c('Spine' = '#f54040', 'Shaft' = '#4da50b')) +
  facet_wrap(~site_type)
  
# site vs int
df.site.type.stat <- df.fret.sites.avg %>%
  group_by(time_interval) %>%
  wilcox_test(int_interval ~ site_type) %>%
  add_significance() %>%
  add_xy_position(fun = "max") 

ggplot() +
  geom_boxplot(data = df.fret.sites.avg,
               aes(x = site_type,
                   y = int_interval,
                   fill = site_type)) +
  geom_point(data = df.fret.sites.avg,
             aes(x = site_type,
                 y = int_interval),
             size=2, shape=21) +
  stat_pvalue_manual(df.site.type.stat, label = 'p.signif',
                     hide.ns = TRUE) +
  theme_classic() +
  theme(legend.position="none") +
  labs(title = 'FRET between site types in different time intervals',
       x = 'Site type',
       y = expression(E[app])) +
  scale_fill_manual(values = c('Spine' = '#f54040', 'Shaft' = '#4da50b')) +
  facet_wrap(~time_interval)


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
ggplot() +
  stat_summary(data = df.hpca.sites,
               aes(x = time, y = int, color=site_type, group = site_type),
               fun = median,
               geom = 'line', size = 0.5) +
  stat_summary(data = df.hpca.sites,
               aes(x = time, y = int, color = site_type, group = site_type),
               fun = median,
               geom = 'point', size = 1) +
  stat_summary(data = df.hpca.sites,
               aes(x = time, y = int,
                   color = site_type, fill = site_type, group = site_type),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median,
               geom = 'ribbon', size = 0, alpha = .25) +
  annotate('text', label = 'NMDA app.', x = 92, y = 0.22, color = 'red') +
  annotate('rect', xmin = 60, xmax = 120, ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = 'red') +
  annotate('text', label = 'Base.', x = 25, y = 0.16, color = 'black') +
  annotate('rect',
           xmin = base.indexes[1] * 10, xmax = rev(base.indexes)[1] * 10,
           ymin = -Inf, ymax = 0.18,
           alpha = 0.075, color = 'black', linewidth = 0) +
  annotate('text', label = 'Mid.', x = 80, y = 0.16, color = 'black') +
  annotate('rect',
           xmin = mid.indexes.h[1] * 10, xmax = rev(mid.indexes.h)[1] * 10,
           ymin = -Inf, ymax = 0.18,
           alpha = 0.075, color = 'black', linewidth = 0) +
  annotate('text', label = 'End', x = 260, y = 0.16, color = 'black') +
  annotate('rect',
           xmin = end.indexes[1] * 10, xmax = rev(end.indexes)[1] * 10,
           ymin = -Inf, ymax = 0.18,
           alpha = 0.075, color = 'black', linewidth = 0) +
  scale_color_manual(values = c('Spine' = '#f54040', 'Shaft' = '#4da50b')) +
  scale_fill_manual(values = c('Spine' = '#f54040', 'Shaft' = '#4da50b')) +
  labs(title = 'HPCA insertions in different ROI types',
       color = 'ROI type',
       fill = 'ROI type',
       x = 'Time, s',
       y = expression(ΔF/F[0])) +
  theme_classic()


# slope
df.hpca.sites.med <- df.hpca.sites %>%
  group_by(site_type, index) %>%
  mutate(int_med = median(int)) %>%
  select(-roi_id, -cell_id, -int) %>%
  distinct() %>%
  ungroup()

o1.h = beast(df.hpca.sites.med$int_med[df.hpca.sites.med$site_type == 'Spine'],
           season = 'none', method='bic')
df.hpca.sites.beast.spine <- data.frame('Time' = seq(0, length(o1.h$trend$slpSgnPosPr)-1) * 10,
                                        'Pos.' = o1.h$trend$slpSgnPosPr,
                                        'Zero' = o1.h$trend$slpSgnZeroPr,
                                        'Neg.' = rep(1, length(o1.h$trend$slpSgnPosPr)) - (o1.h$trend$slpSgnPosPr+o1.h$trend$slpSgnZeroPr)) %>%
  pivot_longer(cols = 'Pos.':'Neg.', names_to = 'SlopeSign', values_to = 'p') %>%
  mutate(SlopeSign = factor(SlopeSign, c('Pos.', 'Zero', 'Neg.'), ordered = TRUE))
ggplot() +
  geom_area(data = df.hpca.sites.beast.spine,
            aes(x = Time, y = p, fill = SlopeSign), alpha = .6) +
  annotate('text', label = 'NMDA app.', x = 90, y = 1.03, color = 'red') +
  annotate('rect', xmin = 60, xmax = 120, ymin = 0, ymax = 1.06,
           alpha = 0.2, fill = 'red') +
  geom_hline(yintercept = 0.59, lty = 2) +
  scale_fill_manual(values = c('Pos.' = '#ff4747',
                               'Zero' = '#49cf36',
                               'Neg.' =  '#3636cf')) +
  scale_x_continuous(limits = c(0,280), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.002,1.06), expand = c(0, 0)) +
  theme_classic() +
  labs(title = 'Probability of slope sign for HPCA insertions, spines sites',
       fill = 'Slope sing',
       x = 'Time, s',
       y = 'Probability')  

o2.h = beast(df.hpca.sites.med$int_med[df.hpca.sites.med$site_type == 'Shaft'],
             season = 'none', method='bic')
df.hpca.sites.beast.shaft <- data.frame('Time' = seq(0, length(o2.h$trend$slpSgnPosPr)-1) * 10,
                                        'Pos.' = o2.h$trend$slpSgnPosPr,
                                        'Zero' = o2.h$trend$slpSgnZeroPr,
                                        'Neg.' = rep(1, length(o2.h$trend$slpSgnPosPr)) - (o2.h$trend$slpSgnPosPr+o2.h$trend$slpSgnZeroPr)) %>%
  pivot_longer(cols = 'Pos.':'Neg.', names_to = 'SlopeSign', values_to = 'p') %>%
  mutate(SlopeSign = factor(SlopeSign, c('Pos.', 'Zero', 'Neg.'), ordered = TRUE))
ggplot() +
  geom_area(data = df.hpca.sites.beast.shaft,
            aes(x = Time, y = p, fill = SlopeSign), alpha = .6) +
  annotate('text', label = 'NMDA app.', x = 90, y = 1.03, color = 'red') +
  annotate('rect', xmin = 60, xmax = 120, ymin = 0, ymax = 1.06,
           alpha = 0.2, fill = 'red') +
  geom_hline(yintercept = 0.0, lty = 2) +
  scale_fill_manual(values = c('Pos.' = '#ff4747',
                               'Zero' = '#49cf36',
                               'Neg.' =  '#3636cf')) +
  scale_x_continuous(limits = c(0,280), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.002,1.06), expand = c(0, 0)) +
  theme_classic() +
  labs(title = 'Probability of slope sign for HPCA insertions, shaft sites',
       fill = 'Slope sing',
       x = 'Time, s',
       y = 'Probability')  


# time intervals
df.hpca.sites.avg <- df.hpca.sites %>%
  select(-time) %>%
  group_by(roi_id, site_type) %>%
  mutate(time_interval = case_when(index %in% base.indexes ~ 'base',
                                   index %in% mid.indexes.h ~ 'mid',
                                   index %in% end.indexes ~ 'end',
                                   .default = 'out')) %>%
  filter(time_interval != 'out') %>%
  droplevels() %>%
  ungroup() %>%
  mutate(time_interval = factor(time_interval, c('base', 'mid', 'end'), ordered = TRUE)) %>%
  group_by(roi_id, site_type, time_interval) %>%
  mutate(int_interval = median(int)) %>%
  select(-index, -int) %>%
  ungroup() %>%
  distinct()

# base vs end
df.base.hpca.zero.stat <- df.hpca.sites.avg %>%
  filter(time_interval == 'end') %>%
  select(-time_interval) %>%
  group_by(site_type) %>%
  wilcox_test(int_interval ~ 1, mu = 0) %>%
  add_significance() %>%
  mutate(y.position = c(0.03,0.02), group2 = c(1,1))

ggplot() +
  geom_boxplot(data = df.hpca.sites.avg %>% filter(time_interval == 'base'),
               aes(x = time_interval,
                   y = int_interval,
                   fill = site_type)) +
  geom_point(data = df.hpca.sites.avg %>% filter(time_interval == 'base'),
             aes(x = time_interval,
                 y = int_interval),
             size=2, shape=21) +
  stat_pvalue_manual(df.base.hpca.zero.stat, label = 'p.signif',
                     hide.ns = TRUE, remove.bracket = TRUE, label.size = 5) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c('Spine' = '#f54040', 'Shaft' = '#4da50b')) +
  scale_y_continuous(limits = c(-0.04, 0.03)) +
  labs(title = 'HPCA amount at the end of observations in different ROI types',
       y = expression(ΔF/F[0])) +
  facet_wrap(~site_type)

# time vs int
df.time.interval.hpca.stat <- df.hpca.sites.avg %>%
  group_by(site_type) %>%
  pairwise_wilcox_test(int_interval ~ time_interval, p.adjust.method = 'BH') %>%
  add_significance() %>%
  add_xy_position(fun = "max") 

ggplot() +
  geom_boxplot(data = df.hpca.sites.avg,
               aes(x = time_interval,
                   y = int_interval,
                   fill = site_type)) +
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
                     hide.ns = TRUE) +
  theme_classic() +
  theme(legend.position="none") +
  labs(title = 'Time courses of HPCA insertions in different ROI types',
       x = 'Time interval',
       y = expression(ΔF/F[0])) +
  scale_fill_manual(values = c('Spine' = '#f54040', 'Shaft' = '#4da50b')) +
  facet_wrap(~site_type)

# site vs int
df.site.type.hpca.stat <- df.hpca.sites.avg %>%
  group_by(time_interval) %>%
  wilcox_test(int_interval ~ site_type) %>%
  add_significance() %>%
  add_xy_position(fun = "max") 

ggplot() +
  geom_boxplot(data = df.hpca.sites.avg,
               aes(x = site_type,
                   y = int_interval,
                   fill = site_type)) +
  geom_point(data = df.hpca.sites.avg,
             aes(x = site_type,
                 y = int_interval),
             size=2, shape=21) +
  stat_pvalue_manual(df.site.type.hpca.stat, label = 'p.signif',
                     hide.ns = TRUE) +
  theme_classic() +
  theme(legend.position="none") +
  labs(title = 'FRET between site types in different time intervals',
       x = 'Site type',
       y = expression(ΔF/F[0])) +
  scale_fill_manual(values = c('Spine' = '#f54040', 'Shaft' = '#4da50b')) +
  facet_wrap(~time_interval)
