# Sasha's data analysis, HPCA+PSD95 FRET control

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

setwd('/home/wisstock/bio_note/projects/PhD/6_2021_NMDAR_plasticity/exp/2026_03_5_HPCA_PSD_FRET_bleach')


##### DF PREPROCESSING #####
df.full <- read.csv('df_combined.csv') %>%
           select(-X) %>%
           mutate_if(is.character, factor) %>%
           mutate(roi = as.factor(roi),
                  roi_id = interaction(id, roi, sep = '_'))
  

df.full.id.summary <- df.full %>%
  filter(base == 'simple', index == 10) %>%
  group_by(id, ch) %>%
  summarise(n_roi = n_distinct(roi_id),
            max = max(dF0_int),
            min = min(dF0_int))

##### PLOT PARAMETERS #####
start_box <- seq(0, 2)
end_box <- seq(117, 120)

ch0.color <- 'green2'
ch1.color <- 'orange2'
ch3.color <- 'red2'
fc.color <- 'blue3'
ea.color <- 'magenta3'


##### EXPLORATION #####
ggplot(data = df.full %>% filter(base == 'simple'),
       aes(x = time, y = abs_int, color = ch, fill = ch,
           linetype = data, shape = data)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  facet_wrap(~id)


ggplot(data = df.full %>% filter(base == 'simple', ch == 'ch0'),
       aes(x = time, y = abs_int, color = data, fill = data,
           group = roi_id)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'line', linewidth = 0.3) +
  # stat_summary(fun = median,
  #              geom = 'point', size = 1) +
  # stat_summary(fun.min = function(z) {quantile(z,0.25)},
  #              fun.max = function(z) {quantile(z,0.75)},
  #              fun = median,
  #              geom = 'ribbon', linewidth = 0, alpha = .15) +
  theme_minimal() +
  facet_wrap(~id)


##### CH PROFILES #####
df.ch <- df.full %>%
  filter(base == 'simple', ch %in% c('ch0', 'ch1', 'ch3'), data == 'target') %>%
  droplevels()

ggplot(data = df.ch,
       aes(x = index, y = dF0_int, colour = ch, fill = ch)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun = median,
               geom = 'line', size = 0.3) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15)


ggplot(data = df.ch,
       aes(x = index, y = abs_int, colour = ch, fill = ch)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun = median,
               geom = 'line', size = 0.3) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15)



df.ch.box <- df.ch %>%
  mutate(box = as.factor(case_when(index %in% start_box ~ 'start',
                         index %in% end_box ~ 'end',
                         .default = 'x'))) %>%
  filter(box != 'x') %>%
  droplevels() %>%
  group_by(ch, box, roi_id) %>%
  mutate(box = factor(box, levels = c('start', 'end'), ordered = TRUE),
         med_int = median(abs_int),
         med_dF = median(dF_int),
         med_dF0 = median(dF0_int)) %>%
  ungroup() %>%
  select(ch, box, med_int, med_dF, med_dF0, roi_id) %>%
  distinct()
  
length(levels(df.ch.box$roi_id))

ggplot(data = df.ch.box,
       aes(x = box, y = med_int, fill = ch)) +
  geom_boxplot() +
  facet_wrap(~ch)


##### Fc NOISE ANALYSIS #####
df.fret <- df.full %>%
  filter(base == 'simple', ch %in% c('Fc', 'Ea')) %>%
  droplevels()


ggplot(data = df.fret %>% filter(ch == 'Fc'),
       aes(x = time, y = abs_int, color = data, fill = data)) +
  geom_hline(yintercept = 0, linetype = 2) +
  stat_summary(fun = median,
               geom = 'point', size = 1) +
  stat_summary(fun = median,
               geom = 'line', size = 0.3) +
  stat_summary(fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               fun = median,
               geom = 'ribbon', linewidth = 0, alpha = .15) +
  ylim(c(-1,10)) +
  theme_minimal()


df.fret.box <- df.fret %>%
  mutate(box = as.factor(case_when(index %in% start_box ~ 'start',
                                   index %in% end_box ~ 'end',
                                   .default = 'x'))) %>%
  filter(box != 'x') %>%
  droplevels() %>%
  group_by(ch, box, roi_id, data) %>%
  mutate(box = factor(box, levels = c('start', 'end'), ordered = TRUE),
         med_int = median(abs_int),
         med_dF = median(dF_int + 0.000001),
         med_dF0 = median(dF0_int + 0.000001)) %>%
  ungroup() %>%
  select(ch, box, med_int, med_dF, med_dF0, roi_id, data) %>%
  distinct()

ggplot(data = df.ch.box %>% filter(ch == 'Fc'),
       aes(x = data, y = med_int)) +
  geom_boxplot()


##### 2D #####
df.wide <- df.full %>%
  filter(base == 'simple', ch %in% c('ch3', 'ch0', 'Fc', 'Ea')) %>%
  select(id, roi, index, abs_int, ch) %>%
  pivot_wider(names_from = ch, values_from = abs_int) %>%
  mutate(rel = ch0/ch3)


ggplot(data = df.wide,
       aes(x = rel, y = Fc, colour = id)) +
  geom_point(alpha = .5)
  scale_x_log10() +
  scale_y_log10()
