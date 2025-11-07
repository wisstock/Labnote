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


df.full <- read.csv('df_processed.csv') %>%
           select(-X) %>%
           mutate_if(is.character, factor) %>%
           mutate(app_time = as.factor(app_time),
                  roi = as.factor(roi),
                  roi_id = interaction(id, roi, sep = '_'),
                  dF_F0_int = dF.F0_int) %>%
           select(-dF.F0_int)

df.full.summary <- df.full %>%
  filter(ch == 'ch0', index == 10, app_time == '60') %>%
  group_by(id, lab_id) %>%
  summarise(n_roi = n_distinct(roi_id),
            max = max(dF_F0_int),
            min = min(dF_F0_int))
  


spines_id <- df.full %>%
  filter(lab_id != 'shaft', lab_id != 'psd_dots') %>%
  droplevels() %>%
  group_by(roi_id) %>%
  summarise(has_all = n_distinct(lab_id) == 2) %>%
  filter(has_all) %>%
  droplevels() %>%
  pull(roi_id)

df.spines <- df.full %>%
  filter(lab_id != 'shaft', lab_id != 'psd_dots') %>%
  filter(roi_id %in% spines_id) %>%
  droplevels() %>%
  group_by(roi_id) %>%
  mutate(dist_um = dist * 0.16) %>%
  droplevels() %>%
  select(-index) %>%
  ungroup()

df.spines.05 <- df.spines %>%
  filter(app_time == '0.5') %>%
  select(-app_time) %>%
  mutate(rel_time = time - 1)

df.spines.60 <- df.spines %>%
  filter(app_time == '60') %>%
  select(-app_time) %>%
  mutate(rel_time = time - 40)


df.shaft.05 <- df.full %>%
  filter(lab_id == 'shaft', app_time == '0.5') %>%
  mutate(dist_um = dist * 0.16,
         rel_time = time - 1) %>%
  droplevels() %>%
  select(-index, -app_time, -lab_id) %>%
  ungroup()

df.shaft.60 <- df.full %>%
  filter(lab_id == 'shaft', app_time == '60') %>%
  mutate(dist_um = dist * 0.16,
         rel_time = time - 40) %>%
  droplevels() %>%
  select(-index, -app_time, -lab_id) %>%
  ungroup()


##### FILTERING #####
low.ids <- df.spines.60 %>%
  filter(base == 'dietrich',
         rel_time == 20,
         (dF_F0_int < 0.05 & lab_id == 'halo' & ch == 'ch0')) %>%
  droplevels()
low.ids.list <- levels(low.ids$roi_id)

low.shafts <- df.shaft.60 %>%
  filter(base == 'simple',
         rel_time == 30,
         (dF_F0_int < 0.05 & ch == 'ch0')) %>%
  droplevels()
low.shafts.list <- levels(low.shafts$roi_id)

remove(low.ids, low.shafts)

bad.cells.list <- c('25_09_11_cell2', '25_09_11_cell1', '25_10_09_cell1')

##### EXPLORATORY ANALYSIS #####
df.spines.summary <- df.spines.60 %>%
  filter(base == 'simple', lab_id == 'psd') %>%
  group_by(id) %>%
  summarise(n_roi = n_distinct(roi_id),
            max = max(abs_int),
            min = min(abs_int))
remove(df.spines.summary)

##### 60 SPINES and SHAFT PDFs #####
ggplot(df.spines.60 %>% filter(! roi_id %in% low.ids.list,
                               base == 'dietrich',
                               rel_time == 10),
       aes(x = dF_F0_int, fill = ch)) +
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_density(alpha = .5) +
  theme_minimal() +
  facet_wrap(~interaction(lab_id,ch), ncol = 2, scales = 'free')

ggplot(df.shaft.60 %>% filter(! roi_id %in% low.ids.list,
                               base == 'dietrich',
                               rel_time == 10),
       aes(x = dF_F0_int, fill = ch)) +
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_density(alpha = .5) +
  theme_minimal() +
  facet_wrap(~ch, ncol = 2, scales = 'free')

##### 0.5 SPINES and SHAFT PDFs #####
ggplot(df.spines.05 %>% filter(! roi_id %in% low.ids.list,
                               base == 'simple',
                               rel_time == 0.25),
       aes(x = dF_F0_int, fill = ch)) +
  geom_vline(xintercept = 0, linetype = 2) + 
  geom_density(alpha = .5) +
  theme_minimal() +
  facet_wrap(~interaction(lab_id,ch), ncol = 2, scales = 'free')

###### 0.5 S PROF #####
ggplot(data = df.spines.05 %>% filter(! roi_id %in% low.ids.list, base == 'simple'),
       aes(x = rel_time, y = dF_int, color = id, fill = id)) +
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

ggplot(data = df.shaft.05 %>% filter(base == 'simple',
                                     !roi_id %in% low.shafts.list),
       aes(x = rel_time, y = dF_int, color = id, fill = id)) +
  geom_vline(xintercept = 0, linetype = 2) +
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


###### 60 S PROF #####
ggplot(data = df.spines.60 %>% filter(base == 'dietrich',
                                      !roi_id %in% low.ids.list,
                                      !id %in% bad.cells.list),
       aes(x = rel_time, y = dF_int, color = id, fill = id)) +
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

ggplot(data = df.shaft.60 %>% filter(base == 'dietrich',
                                     !roi_id %in% low.shafts.list,
                                     !id %in% bad.cells.list),
       aes(x = rel_time, y = abs_int, color = id, fill = id)) +
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


###### 60 S HALO vs PSD #####
df.spines.60.wide <- df.spines.60 %>%
  filter(base == 'simple', ch %in% c('Fc', 'ch0')) %>%
  droplevels() %>%
  mutate(lab_ch = interaction(lab_id, ch)) %>%
  select(lab_ch, rel_time, roi_id, id, dF_F0_int) %>%
  distinct() %>%
  pivot_wider(names_from = lab_ch, values_from = dF_F0_int)

df.avg.roi <- df.by.roi %>%
  filter(app_factor %in% c('20'), dist_group != 'I') %>%  # 
  group_by(rel_time) %>%
  summarise(psd_mean = mean(psd, na.rm = TRUE),
            oreol_mean = mean(oreol, na.rm = TRUE),
            psd_se = sd(psd, na.rm = TRUE) / sqrt(n()),
            oreol_se = sd(oreol, na.rm = TRUE) / sqrt(n()),
            n = n(),
            .groups = 'drop')


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


ggplot(data = df.spines.60.wide %>% filter(rel_time == 20, !roi_id %in% low.ids.list),
       aes(x = psd.ch0, y = halo.Fc, color = roi_id)) +
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
                                      base == 'simple'),
       aes(x = rel_time, y = abs_int, color = lab_id, fill = lab_id)) +
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
                                      base == 'dietrich'),
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

