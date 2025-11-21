# Sasha's data analysis, HPCA+AP2 FRET

library(dplyr)
library(tidyr)
library(purrr)
library(rstatix)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(cowplot)
library(ggsci)
library(introdataviz)

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2025_11_03_AP2_FRET_Bao')


##### DF PREPROCESSING #####
# E_A
epsilon.A.405 <- 1.3
epsilon.D.405 <- 94.3
a <- 0.0135
e.div.e <- epsilon.A.405/epsilon.D.405


df.full <- read.csv('df_combined.csv') %>%
           select(-X) %>%
           mutate_if(is.character, factor) %>%
           mutate(roi = as.factor(roi),
                  roi_id = interaction(id, roi, sep = '_')) %>%
           group_by(id, roi, index)

df.full.id.summary <- df.full %>%
  filter(ch == 'ch0', index == 10) %>%
  group_by(id, lab_id) %>%
  summarise(n_roi = n_distinct(roi_id),
            max = max(dF0_int),
            min = min(dF0_int))

df.full.lab.summary <- df.full %>%
  filter(ch == 'ch0', index == 10) %>%
  group_by(lab_id) %>%
  summarise(n_roi = n_distinct(roi_id),
            max = max(dF0_int),
            min = min(dF0_int))
  
# SPINES FILTERING
spines_id <- df.full %>%
  filter(lab_id %in% c('psd', 'halo')) %>%
  droplevels() %>%
  group_by(roi_id) %>%
  summarise(has_all = n_distinct(lab_id) == 2) %>%
  filter(has_all) %>%
  droplevels() %>%
  pull(roi_id)

df.spines <- df.full %>%
  filter(lab_id %in% c('psd', 'halo')) %>%
  filter(roi_id %in% spines_id) %>%
  droplevels() %>%
  group_by(roi_id) %>%
  mutate(dist_um = dist * 0.16) %>%
  droplevels() %>%
  select(-index) %>%
  ungroup()


df.shaft <- df.full %>%
  filter(lab_id == 'shaft') %>%
  mutate(dist_um = dist * 0.16,
         rel_time = time - 40) %>%
  droplevels() %>%
  select(-index, -lab_id) %>%
  ungroup()


df.ccp <- df.full %>%
  filter(lab_id == 'ccp') %>%
  mutate(dist_um = dist * 0.16,
         rel_time = time - 40) %>%
  droplevels() %>%
  select(-index, -lab_id) %>%
  ungroup()


##### FILTERING #####
low.spines <- df.spines %>%
  filter(base == 'dietrich',
         rel_time == 20,
         (dF0_int < 0.075 & lab_id == 'halo' & ch == 'ch0')) %>%
  droplevels() %>%
  pull(roi_id)

low.shafts <- df.shaft %>%
  filter(base == 'simple',
         rel_time == 20,
         (dF0_int < 0.075 & ch == 'ch0')) %>%
  droplevels() %>%
  pull(roi_id)

# low.ccps <- df.ccp %>%
#   filter(base == 'simple',
#          rel_time == 20,
#          (dF0_int > 0.075 & ch == 'ch0')) %>%
#   droplevels() %>%
#   pull(roi_id)

bad.cells.list <- c('25_09_10_cell1')


##### SPINES and SHAFT PDFs #####
ggplot(df.spines %>% filter(! roi_id %in% low.spines,
                            base == 'dietrich',
                            rel_time == 20),
       aes(x = dF0_int, fill = ch)) +
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_density(alpha = .5) +
  theme_minimal() +
  facet_wrap(~interaction(lab_id,ch), ncol = 2, scales = 'free')

ggplot(df.shaft %>% filter(! roi_id %in% low.shafts,
                           base == 'dietrich',
                           rel_time == 30),
       aes(x = dF0_int, fill = ch)) +
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_density(alpha = .5) +
  theme_minimal() +
  facet_wrap(~ch, ncol = 2, scales = 'free')


###### PROF PER CELL #####
# spines
ggplot(data = df.spines %>% filter(base == 'dietrich',
                                   !roi_id %in% low.spines,
                                   !id %in% bad.cells.list),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .05) +
  theme_minimal() +
  facet_wrap(~interaction(lab_id,ch), ncol = 2, scales = 'free_y')

# shaft
ggplot(data = df.shaft %>% filter(base == 'dietrich',
                                  !roi_id %in% low.shafts,
                                  !id %in% bad.cells.list),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .05) +
  theme_minimal() +
  facet_wrap(~ch, ncol = 2, scales = 'free_y')

# CCP
ggplot(data = df.ccp %>% filter(base == 'dietrich',
                                  !id %in% bad.cells.list),
       aes(x = rel_time, y = dF0_int, color = id, fill = id)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .05) +
  theme_minimal() +
  facet_wrap(~ch, ncol = 2, scales = 'free_y')

##### PROF AVG #####
ggplot(data = df.spines %>% filter(base == 'dietrich',
                                   !roi_id %in% low.spines,
                                   !id %in% bad.cells.list),
       aes(x = rel_time, y = dF0_int, color = lab_id, fill = lab_id)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .05) +
  theme_minimal() +
  facet_wrap(~ch, ncol = 2, scales = 'free_y')

ggplot(data = df.shaft %>% filter(base == 'dietrich',
                                  !roi_id %in% low.shafts,
                                  !id %in% bad.cells.list),
       aes(x = rel_time, y = dF0_int, color = ch, fill = ch)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 60, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .05) +
  theme_minimal()


###### HALO vs PSD #####
df.spines.wide.f <- df.spines %>%
  filter(base == 'dietrich', ch %in% c('ch0', 'ch3', 'Fc')) %>%
  droplevels() %>%
  mutate(lab_ch = interaction(lab_id, ch)) %>%
  select(lab_ch, rel_time, roi_id, id, dF0_int) %>%
  distinct() %>%
  pivot_wider(names_from = lab_ch, values_from = dF0_int)

df.spines.wide.e <- df.spines %>%
  filter(base == 'dietrich', ch %in% c('Eapp')) %>%
  droplevels() %>%
  mutate(lab_ch = interaction(lab_id, ch)) %>%
  select(lab_ch, rel_time, roi_id, id, abs_int) %>%
  distinct() %>%
  pivot_wider(names_from = lab_ch, values_from = abs_int)

df.spines.wide <- cbind(df.spines.wide.f, df.spines.wide.e[, 4:5])
remove(df.spines.wide.f, df.spines.wide.e)

df.avg.spines <- df.spines.wide %>%
  group_by(rel_time) %>%
  summarise(psd_i_mean = mean(psd.ch0, na.rm = TRUE),
            psd_i_se = sd(psd.ch0, na.rm = TRUE) / sqrt(n()),
            halo_i_mean = mean(halo.ch0, na.rm = TRUE),
            halo_i_se = sd(halo.ch0, na.rm = TRUE) / sqrt(n()),
            psd_e_mean = mean(psd.Fc, na.rm = TRUE),
            psd_e_se = sd(psd.Fc, na.rm = TRUE) / sqrt(n()),
            halo_e_mean = mean(halo.Fc[is.finite(halo.Fc)], na.rm = TRUE),
            halo_e_se = sd(halo.Fc[is.finite(halo.Fc)], na.rm = TRUE) / sqrt(n()),
            n = n(),
            .groups = 'drop')

track.time <- seq(0,120)

ggplot(data = df.avg.spines %>% filter(rel_time %in% track.time),
       aes(x = psd_i_mean, y = halo_i_mean, color = rel_time)) +
  # geom_vline(xintercept = 0, linetype = 2) +
  # geom_hline(yintercept = 0, linetype = 2) +
  geom_path() +
  geom_point()
  geom_linerange(aes(xmin = psd_i_mean - psd_i_se,
                     xmax = psd_i_mean + psd_i_se),
                 size = 0.3) +
  geom_linerange(aes(ymin = halo_e_mean - halo_e_se,
                     ymax = halo_e_mean + halo_e_se,
                     color = rel_time),
                 size = 0.3)

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


ggplot(data = df.spines.60.wide %>% filter(rel_time == 30, !roi_id %in% low.ids.list),
       aes(x = psd.ch0, y = halo.ch0, color = roi_id)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = .3) +
  theme(legend.position = "none") +
  facet_wrap(~id)

ggplot(data = df.spines.60.wide %>% filter(rel_time >= 0, rel_time <= 100),
       aes(x = psd.ch0, y = halo.ch0, color = rel_time)) +
  geom_point(alpha = .3) +
  facet_wrap(~id)


#### CELL SUMMARY PLOT #####
ggplot(data = df.spines.60 %>% filter(!roi_id %in% low.ids.list,
                                      ch %in% c('Eapp', 'Fc'),  # c('ch0', 'ch3'),
                                      base == 'simple',
                                      !id %in% bad.cells.list),
       aes(x = rel_time, y = dF_F0_int, color = lab_id, fill = lab_id)) +
  geom_vline(xintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .1) +
  theme_minimal() +
  facet_wrap(~ch, ncol = 2, scales = 'free_y')

ggplot(data = df.spines.60 %>% filter(!roi_id %in% low.ids.list,
                                      ch %in% c('ch0', 'ch3'),
                                      base == 'dietrich',
                                      !id %in% bad.cells.list),
       aes(x = rel_time, y = dF_F0_int, color = lab_id, fill = lab_id)) +
  geom_vline(xintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .1) +
  theme_minimal() +
  facet_wrap(~ch, ncol = 2)

# shaft
ggplot(data = df.shaft.60 %>% filter(!roi_id %in% low.shafts.list,
                                      ch %in% c('ch0', 'ch3'),
                                      base == 'dietrich',
                                     !id %in% c('25_09_11_cell2', '25_09_11_cell1')),
       aes(x = rel_time, y = dF_F0_int, color = ch, fill = ch)) +
  geom_vline(xintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .1) +
  theme_minimal() +
  facet_wrap(~ch, ncol = 2)

ggplot(data = df.shaft.60 %>% filter(!roi_id %in% low.shafts.list,
                                     ch %in% c('Fc', 'Eapp'),
                                     base == 'simple',
                                     !id %in% c('25_09_11_cell2', '25_09_11_cell1')),
       aes(x = rel_time, y = dF_int, color = ch, fill = ch)) +
  geom_vline(xintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .1) +
  theme_minimal() +
  facet_wrap(~ch, ncol = 2, scale = 'free')

